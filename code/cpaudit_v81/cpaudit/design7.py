"""Phase 6 Task 6.4 design (v7): every split, fold and partition, built once at `setup7`,
stored on Drive and hashed into the pre-registration. Fitting steps read them back; nothing here
uses a target value except the pre-specified scope rule (band gaps > 0 = non-metals).

Regime B (unseen space group), MatFold criterion (Witman & Schindler 2025), nested (K = 10, L = 5):
  * eligible groups = MatFold space groups with >= min_group_size rows in the target's scope;
    rarer groups are never held out and always train (MatFold fraction_lower_limit rule);
  * eligible groups are split into two halves (A, B) for nested tuning; each half is split into
    5 outer folds; repeat 2 re-draws the outer folds inside the same halves;
  * inner folds: 5 folds by space group over the outer training set's eligible groups (MatFold
    "inner = outer criterion"); rare groups train every inner model and calibrate none;
  * random-inner contrast: the same calibration rows, folds random by reduced formula;
  * polymorph exclusion: training drops every row whose reduced formula occurs in the held-out
    groups of the fold;
  * LOCSO: one crystal system per outer fold, inner leave-one-crystal-system-out.
Regime A (in distribution): 5 outer folds x 2 repeats grouped by reduced formula and stratified by
  space group; CV+ folds and a 75/25 proper/calibration split inside each training pool.
H3: 10 partitions; the calibration and test groups come from one half, the other half trains.
"""
from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import pandas as pd

TARGETS = {
    "formation": dict(col="formation_energy_peratom", scope="all", label="formation energy per atom (eV/atom)"),
    "bandgap": dict(col="optb88vdw_bandgap", scope="positive", label="OptB88vdW band gap, non-metals (eV)"),
    "tbmbj": dict(col="mbj_bandgap", scope="positive", label="TBmBJ band gap, non-metals (eV)"),
}


def scope_mask(F: pd.DataFrame, target: str) -> np.ndarray:
    y = pd.to_numeric(F[TARGETS[target]["col"]], errors="coerce").to_numpy(float)
    m = np.isfinite(y)
    if TARGETS[target]["scope"] == "positive":
        m &= y > 0
    return m


def balanced_assign(items, sizes, n_bins: int, rng) -> np.ndarray:
    """Items in random order; each goes to the bin with the fewest items, then the fewest rows."""
    items = np.asarray(items)
    sizes = np.asarray(sizes, float)
    order = rng.permutation(len(items))
    cnt, tot = np.zeros(n_bins, int), np.zeros(n_bins)
    out = np.empty(len(items), int)
    for i in order:
        b = int(np.lexsort((np.arange(n_bins), tot, cnt))[0])
        out[i] = b
        cnt[b] += 1
        tot[b] += sizes[i]
    return out


def _seed(base, *parts):
    h = hashlib.sha256(json.dumps([int(base)] + [str(p) for p in parts]).encode()).hexdigest()
    return int(h[:8], 16)


def build_target_design(F: pd.DataFrame, L: pd.DataFrame, target: str, m7: dict) -> dict:
    """All mappings for one target. Row sets are derived from these by `rowsets` below."""
    base = int(m7["seed"])
    min_size = int(m7["min_group_size"])
    mask = scope_mask(F, target)
    rows = np.flatnonzero(mask)
    grp = L.sga01_number.to_numpy().astype(int)[rows]
    crys = np.asarray(L.sga01_crys.astype(str).str.lower().to_numpy()[rows], dtype=str)
    form = np.asarray(F.reduced_formula.astype(str).to_numpy()[rows], dtype=str)
    ug, cnt = np.unique(grp, return_counts=True)
    elig = ug[cnt >= min_size]
    size = dict(zip(ug, cnt))
    D = dict(target=target, rows=rows, grp=grp, crys=crys, form=form, eligible=elig,
             elig_size=np.array([size[g] for g in elig]))
    rng = np.random.default_rng(_seed(base, target, "halves"))
    half = balanced_assign(elig, D["elig_size"], 2, rng)
    D["half"] = half
    n_outer_half = int(m7["outer_folds_per_half"])
    for rep in (1, 2):
        outer = np.empty(len(elig), int)
        for h in (0, 1):
            idx = np.flatnonzero(half == h)
            r = np.random.default_rng(_seed(base, target, "outer", rep, h))
            outer[idx] = h * n_outer_half + balanced_assign(elig[idx], D["elig_size"][idx], min(n_outer_half, len(idx)), r)
        D[f"outer_r{rep}"] = outer
    # inner group folds per (rep, outer fold): over the eligible groups that train
    L_in = int(m7["inner_folds"])
    for rep in (1, 2):
        inner = np.full((2 * n_outer_half, len(elig)), -1, int)
        for o in range(2 * n_outer_half):
            tr = np.flatnonzero(D[f"outer_r{rep}"] != o)
            r = np.random.default_rng(_seed(base, target, "inner", rep, o))
            inner[o, tr] = balanced_assign(elig[tr], D["elig_size"][tr], L_in, r)
        D[f"inner_r{rep}"] = inner
    # random-inner folds: one fold per reduced formula (seeded), per outer fold of repeat 1
    uf = np.unique(form)
    D["formulas"] = uf
    rf = np.empty((2 * n_outer_half, len(uf)), np.int8)
    for o in range(2 * n_outer_half):
        r = np.random.default_rng(_seed(base, target, "randinner", o))
        rf[o] = r.permutation(np.arange(len(uf)) % L_in)
    D["randinner_formula_fold"] = rf
    # tuning folds: for half h, 3 group folds over the other half's eligible groups
    for h in (0, 1):
        other = np.flatnonzero(half != h)
        r = np.random.default_rng(_seed(base, target, "tune", h))
        D[f"tune_fold_h{h}"] = balanced_assign(elig[other], D["elig_size"][other], int(m7["tuning"]["inner_folds"]), r)
    # crystal systems (LOCSO)
    D["crystal_systems"] = np.array(sorted(set(crys)))
    # Regime A (in distribution): stratified, formula-grouped folds; stored as arrays
    from sklearn.model_selection import StratifiedGroupKFold
    strat = np.where(np.isin(grp, elig), grp, -1)
    fcode = np.searchsorted(uf, form)
    A_outer = np.full((2, len(rows)), -1, np.int8)
    A_cv = np.full((10, len(rows)), -1, np.int8)
    A_proper = np.zeros((10, len(rows)), bool)
    for rep in range(2):
        sg = StratifiedGroupKFold(5, shuffle=True, random_state=_seed(base, target, "A_outer", rep) % (2 ** 31))
        for k, (_, te) in enumerate(sg.split(np.zeros(len(rows)), strat, fcode)):
            A_outer[rep, te] = k
        for k in range(5):
            s = rep * 5 + k
            pool = np.flatnonzero(A_outer[rep] != k)
            sg2 = StratifiedGroupKFold(5, shuffle=True, random_state=_seed(base, target, "A_cv", s) % (2 ** 31))
            for j, (_, v) in enumerate(sg2.split(np.zeros(len(pool)), strat[pool], fcode[pool])):
                A_cv[s, pool[v]] = j
            sg3 = StratifiedGroupKFold(4, shuffle=True, random_state=_seed(base, target, "A_cal", s) % (2 ** 31))
            _, cal = next(sg3.split(np.zeros(len(pool)), strat[pool], fcode[pool]))
            prop = np.ones(len(pool), bool)
            prop[cal] = False
            A_proper[s, pool] = prop
    D.update(A_outer=A_outer, A_cv=A_cv, A_proper=A_proper)
    # Regime A tuning folds: 3 formula-grouped folds over split 0's proper rows
    r = np.random.default_rng(_seed(base, target, "A_tune"))
    D["A_tune_formula_fold"] = r.permutation(np.arange(len(uf)) % int(m7["tuning"]["inner_folds"]))
    return D


# ----------------------------------------------------------------------------- derived row sets
def B_rowsets(D: dict, rep: int, o: int, design: str) -> dict:
    """Positions (into D['rows']) for one Regime B unit, and the inner fold of each calibration row."""
    grp, elig = D["grp"], D["eligible"]
    outer = D[f"outer_r{rep}"]
    test_groups = elig[outer == o]
    is_test = np.isin(grp, test_groups)
    is_elig = np.isin(grp, elig)
    train = np.flatnonzero(~is_test)
    if design == "polyex":
        bad = np.isin(D["form"], np.unique(D["form"][is_test]))
        train = np.flatnonzero(~is_test & ~bad)
    g2f = dict(zip(elig, D[f"inner_r{rep}"][o]))
    cal = train[is_elig[train]]
    if design == "randinner":
        fcode = np.searchsorted(D["formulas"], D["form"][cal])
        cal_fold = D["randinner_formula_fold"][o][fcode].astype(int)
    else:
        cal_fold = np.array([g2f[g] for g in grp[cal]], int)
    rare = train[~is_elig[train]]
    return dict(test=np.flatnonzero(is_test), train=train, cal=cal, cal_fold=cal_fold, rare=rare,
                test_groups=test_groups, n_folds=int(D[f"inner_r{rep}"].max()) + 1,
                half=int(D["half"][outer == o][0]))


def LOCSO_rowsets(D: dict, o: int) -> dict:
    cs = D["crystal_systems"]
    c = cs[o]
    is_test = D["crys"] == c
    train = np.flatnonzero(~is_test)
    others = [x for x in cs if x != c]
    cal_fold = np.searchsorted(np.array(others), D["crys"][train])
    return dict(test=np.flatnonzero(is_test), train=train, cal=train, cal_fold=cal_fold.astype(int),
                rare=np.array([], int), test_groups=np.array([c]), n_folds=len(others), half=0)


def A_rowsets(D: dict, s: int) -> dict:
    rep, k = divmod(s, 5)
    test = np.flatnonzero(D["A_outer"][rep] == k)
    pool = np.flatnonzero(D["A_outer"][rep] != k)
    proper = pool[D["A_proper"][s][pool]]
    cal = pool[~D["A_proper"][s][pool]]
    return dict(test=test, pool=pool, cv_fold=D["A_cv"][s][pool].astype(int), proper=proper, cal=cal)


def H3_partition(D: dict, p: int, m7: dict) -> dict:
    h = p % 2
    pool_groups = D["eligible"][D["half"] == h]
    is_pool = np.isin(D["grp"], pool_groups)
    return dict(half=h, pool_groups=pool_groups, pool=np.flatnonzero(is_pool), train=np.flatnonzero(~is_pool),
                n_cal=int(m7["h3"]["n_cal_groups"][D["target"]]),
                seed=_seed(int(m7["seed"]), D["target"], "h3", p))


def tune_rowsets(D: dict, which: str) -> tuple[np.ndarray, list]:
    """Rows and (train, valid) folds (positions into the returned rows) for one tuning run."""
    grp, elig = D["grp"], D["eligible"]
    if which in ("h0", "h1"):
        h = int(which[1])
        excl = elig[D["half"] == h]
        rows = np.flatnonzero(~np.isin(grp, excl))
        other = elig[D["half"] != h]
        g2f = dict(zip(other, D[f"tune_fold_{which}"]))
        f = np.array([g2f.get(g, -1) for g in grp[rows]])
    else:   # Regime A: split 0's proper rows, formula-grouped
        rows = A_rowsets(D, 0)["proper"]
        fcode = np.searchsorted(D["formulas"], D["form"][rows])
        f = D["A_tune_formula_fold"][fcode]
    folds = [(np.flatnonzero(f != k), np.flatnonzero(f == k)) for k in range(int(f.max()) + 1)]
    return rows, folds


# ----------------------------------------------------------------------------- checks
def check_design(D: dict, m7: dict) -> dict:
    """Invariants of the design; raises on any violation. Returns a summary."""
    grp, elig = D["grp"], D["eligible"]
    n_out = 2 * int(m7["outer_folds_per_half"])
    out = {}
    for rep in (1, 2):
        outer = D[f"outer_r{rep}"]
        assert set(np.unique(outer)) <= set(range(n_out))
        for o in range(n_out):
            gs = elig[outer == o]
            if len(gs) == 0:
                continue
            assert len(set(D["half"][outer == o])) == 1, "an outer fold mixes halves"
            for design in ("primary", "randinner", "polyex"):
                R = B_rowsets(D, rep, o, design)
                assert not np.isin(grp[R["train"]], gs).any(), "held-out group in training"
                assert np.isin(grp[R["test"]], gs).all()
                assert set(np.unique(grp[R["rare"]])).isdisjoint(set(elig)), "eligible group among rare rows"
                assert not np.isin(grp[R["cal"]], gs).any()
                if design == "polyex":
                    assert not np.isin(D["form"][R["train"]], D["form"][R["test"]]).any(), "polymorph left in training"
                if design != "randinner":
                    for k in range(R["n_folds"]):
                        assert len(set(grp[R["cal"][R["cal_fold"] == k]]) & set(grp[R["cal"][R["cal_fold"] != k]])) == 0
        # every eligible group held out exactly once per repeat
        assert len(np.unique(outer)) == min(n_out, len(elig))
    # tuning never sees the half it is used for
    for h in (0, 1):
        rows, folds = tune_rowsets(D, f"h{h}")
        assert not np.isin(grp[rows], elig[D["half"] == h]).any(), "tuning saw its own half"
    # Regime A: formula-disjoint test / pool, proper / cal
    for s in range(10):
        R = A_rowsets(D, s)
        assert set(D["form"][R["test"]]).isdisjoint(set(D["form"][R["pool"]])), "formula in test and pool"
        assert set(D["form"][R["proper"]]).isdisjoint(set(D["form"][R["cal"]])), "formula in proper and cal"
        assert (D["A_cv"][s][R["pool"]] >= 0).all()
    # H3: pool groups from one half; training excludes them
    for p in range(int(m7["h3"]["n_partitions"]) if D["target"] in m7["h3"]["n_cal_groups"] else 0):
        P = H3_partition(D, p, m7)
        assert not np.isin(grp[P["train"]], P["pool_groups"]).any()
        assert len(P["pool_groups"]) > P["n_cal"], "too few pool groups for the calibration size"
    out["n_rows"] = int(len(D["rows"]))
    out["n_eligible"] = int(len(elig))
    out["rows_in_eligible"] = int(np.isin(grp, elig).sum())
    out["half_sizes_groups"] = [int((D["half"] == h).sum()) for h in (0, 1)]
    out["half_sizes_rows"] = [int(D["elig_size"][D["half"] == h].sum()) for h in (0, 1)]
    out["outer_groups_per_fold_r1"] = [int((D["outer_r1"] == o).sum()) for o in range(n_out)]
    out["outer_rows_per_fold_r1"] = [int(D["elig_size"][D["outer_r1"] == o].sum()) for o in range(n_out)]
    out["crystal_systems"] = [str(c) for c in D["crystal_systems"]]
    return out


def content_hash(D: dict) -> str:
    h = hashlib.sha256()
    for k in sorted(D):
        v = D[k]
        h.update(k.encode())
        if isinstance(v, np.ndarray):
            h.update(str(v.dtype).encode())
            h.update(np.ascontiguousarray(v).tobytes() if v.dtype.kind != "U" else "\x00".join(v.ravel().tolist()).encode())
        else:
            h.update(str(v).encode())
    return h.hexdigest()


def save_design(path: str, D: dict):
    from .config import atomic_save_npz
    arrays = {k: (v if isinstance(v, np.ndarray) else np.array(v)) for k, v in D.items()}
    atomic_save_npz(path, **arrays)


def load_design(path: str) -> dict:
    with np.load(path, allow_pickle=False) as Z:
        D = {k: Z[k] for k in Z.files}
    D["target"] = str(D["target"])
    return D
