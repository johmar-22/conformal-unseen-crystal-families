"""Phase 6, Task 6.1: JARVIS-DFT data audit (counts only; decides the JARVIS scope before any model).

Plan (Section 6, Task 6.1): load dft_3d; fields atoms, spg_number, formation_energy_peratom,
optb88vdw_bandgap, mbj_bandgap ("na" -> NaN); filter to oxides (O present and the only anion);
remove duplicates by (reduced formula, spg_number) keeping the lowest-energy entry; count
structures per space group and keep groups with >= 200; if fewer than about 20 space groups
survive, use all chemistries; count the TBmBJ subset.

What this module does, in order:
 1. Downloads the pinned figshare release(s) (the files jarvis-tools' `data("dft_3d")` loads) to
    Drive (phase06_v1/inputs/) once, and records size and SHA-256.
 2. Builds a slim entry table (no structures): reduced formula from the atoms' element list,
    element set, space group, crystal system, presence of each target. Target VALUES are not
    summarised anywhere; only formation energy is used, to keep the lowest-energy duplicate.
 3. Classifies oxides by two rules (config: primary = oxidation state; sensitivity = element rule).
 4. De-duplicates within each scope by (reduced formula, spg_number), keeping the lowest formation
    energy per atom (for one composition this equals the lowest total energy per atom), ties by
    total energy and then the JARVIS id.
 5. Counts per space group (and crystal system): entries, entries with a TBmBJ gap, distinct
    chemical systems, and the share of entries whose formula also occurs in another space group
    (relevant to the planned polymorph-exclusion sensitivity check, Task 6.3).
 6. Applies the decision rule written in the config before the audit is run, and writes the
    eligible space groups of the chosen scope for Phase 6.
"""
from __future__ import annotations

import collections
import hashlib
import json
import os
import time
import warnings
import zipfile

import numpy as np
import pandas as pd

from .config import (atomic_write_bytes, atomic_write_csv, atomic_write_json, atomic_write_text, code_hash,
                     config_hash, get_paths, sha256_file)
from .phase5 import _gz

FIELDS_NUM = ("formation_energy_peratom", "optb88vdw_bandgap", "mbj_bandgap", "optb88vdw_total_energy", "ehull")
CRYSTAL_SYSTEMS = ["triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"]


# ----------------------------------------------------------------------------- download / load
def fetch_release(rel: dict, inputs_dir: str) -> dict:
    """Download <file>.zip from the pinned URL to Drive once; return path, size and SHA-256."""
    import requests
    os.makedirs(inputs_dir, exist_ok=True)
    path = os.path.join(inputs_dir, rel["file"] + ".zip")
    if not os.path.exists(path):
        print(f"  downloading {rel['name']} ({rel['url']}) -> {path}", flush=True)
        part = path + ".part"
        with requests.get(rel["url"], stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(part, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        with zipfile.ZipFile(part) as z:
            if rel["file"] not in z.namelist():
                raise RuntimeError(f"{rel['url']} does not contain {rel['file']} (found {z.namelist()})")
        os.replace(part, path)
    sha = sha256_file(path)
    if rel.get("sha256") and rel["sha256"] != sha:
        raise RuntimeError(f"{path}: SHA-256 {sha} differs from the pinned {rel['sha256']}")
    return dict(name=rel["name"], file=rel["file"], url=rel["url"], path=path, bytes=os.path.getsize(path), sha256=sha)


def load_release(info: dict) -> list[dict]:
    with zipfile.ZipFile(info["path"]) as z:
        return json.loads(z.read(info["file"]))


def _num(v):
    if v is None or (isinstance(v, str) and v.strip().lower() in ("na", "nan", "none", "")):
        return np.nan
    try:
        x = float(v)
    except (TypeError, ValueError):
        return np.nan
    return x if np.isfinite(x) else np.nan


def entry_table(entries: list[dict]) -> pd.DataFrame:
    from pymatgen.core import Composition
    rows = []
    for e in entries:
        els = e["atoms"]["elements"]
        comp = Composition(dict(collections.Counter(els)))
        rf = comp.reduced_formula
        syms = sorted({el.symbol for el in comp.elements})
        r = dict(jid=e["jid"], formula_field=e.get("formula"), reduced_formula=rf, chemsys="-".join(syms),
                 n_elements=len(syms), nat=len(els), spg_number=_num(e.get("spg_number")),
                 spg_symbol=e.get("spg_symbol"), crys=e.get("crys"), dimensionality=e.get("dimensionality"),
                 func=e.get("func"), has_O=("O" in syms and len(syms) > 1))
        for k in FIELDS_NUM:
            r[k] = _num(e.get(k))
        rows.append(r)
    T = pd.DataFrame(rows)
    bad = T.spg_number.isna()
    if bad.any():
        print(f"  {int(bad.sum())} entries without a space-group number are excluded", flush=True)
    T = T[~bad].copy()
    T["spg_number"] = T.spg_number.astype(int)
    T["jid_num"] = pd.to_numeric(T.jid.str.extract(r"(\d+)$")[0], errors="coerce")
    return T


# ----------------------------------------------------------------------------- oxide rules
def _negatives(formula: str):
    """Elements with a negative oxidation state in the most probable charge-balanced assignment
    (pymatgen oxi_state_guesses, common oxidation states, ICSD-based probabilities). None if no
    charge-balanced assignment exists. Large formulas are reduced with a growing max_sites."""
    from pymatgen.core import Composition
    c = Composition(formula)
    for ms in (-50, -100, -200, None):
        try:
            g = c.oxi_state_guesses(max_sites=ms) if ms is not None else c.oxi_state_guesses()
            break
        except ValueError:
            continue
    else:
        return None
    if not g:
        return None
    return tuple(sorted(el for el, v in g[0].items() if v < 0))


def classify_oxides(T: pd.DataFrame, other_anions) -> pd.DataFrame:
    T = T.copy()
    other = set(other_anions)
    T["oxide_element_rule"] = T.has_O & T.chemsys.map(lambda s: not (set(s.split("-")) - {"O"}) & other)
    forms = T.loc[T.has_O, "reduced_formula"].unique()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        neg = {f: _negatives(f) for f in forms}
    status = []
    for has_o, f in zip(T.has_O, T.reduced_formula):
        if not has_o:
            status.append("no_O")
        elif neg[f] is None:
            status.append("no_charge_balanced_assignment")
        elif neg[f] == ("O",):
            status.append("oxide")
        else:
            status.append("other_anion:" + ",".join(neg[f]))
    T["oxide_status"] = status
    T["oxide_oxidation_rule"] = T.oxide_status == "oxide"
    return T


# ----------------------------------------------------------------------------- dedup / counts
def dedup(S: pd.DataFrame) -> pd.DataFrame:
    """Keep the lowest formation energy per atom per (reduced formula, space group); ties broken by
    total energy per atom, then the JARVIS id (NaN energies sort last)."""
    S = S.sort_values(["formation_energy_peratom", "optb88vdw_total_energy", "jid_num"], na_position="last",
                      kind="mergesort")
    return S.drop_duplicates(["reduced_formula", "spg_number"], keep="first")


def group_counts(S: pd.DataFrame, min_n: int) -> pd.DataFrame:
    """Per space group of a de-duplicated scope."""
    n_spg_of_formula = S.groupby("reduced_formula").spg_number.nunique()
    poly = S.reduced_formula.map(n_spg_of_formula) > 1
    S = S.assign(_poly=poly.to_numpy(), _mbj=S.mbj_bandgap.notna().to_numpy())
    G = S.groupby("spg_number").agg(n=("jid", "size"), spg_symbol=("spg_symbol", "first"), crys=("crys", "first"),
                                    n_mbj=("_mbj", "sum"), n_chemsys=("chemsys", "nunique"),
                                    share_formula_in_other_spg=("_poly", "mean")).reset_index()
    top = S.groupby(["spg_number", "chemsys"]).size().groupby("spg_number").max()
    G["largest_chemsys_share"] = G.spg_number.map(top) / G.n
    G["eligible"] = G.n >= min_n
    return G.sort_values("n", ascending=False).reset_index(drop=True)


def scopes(T: pd.DataFrame) -> dict:
    return {"all_chemistries": T, "oxides_oxidation_rule": T[T.oxide_oxidation_rule],
            "oxides_element_rule": T[T.oxide_element_rule]}


def audit_release(T: pd.DataFrame, a: dict, release: str):
    min_n = int(a["min_group_size"])
    sp_rows, th_rows, cs_rows, sc_rows = [], [], [], []
    eligible = {}
    for scope, S0 in scopes(T).items():
        S = dedup(S0)
        G = group_counts(S, min_n)
        G.insert(0, "scope", scope)
        G.insert(0, "release", release)
        sp_rows.append(G)
        # TBmBJ secondary target: de-duplicate within the entries that have a TBmBJ gap
        Sm = dedup(S0[S0.mbj_bandgap.notna()])
        nm = Sm.groupby("spg_number").size()
        K = int(G.eligible.sum())
        eligible[scope] = G[G.eligible]
        sc_rows.append(dict(release=release, scope=scope, n_raw=len(S0), n_dedup=len(S), n_removed=len(S0) - len(S),
                            n_spacegroups=int(G.shape[0]), K_eligible=K,
                            share_of_entries_in_eligible=float(G[G.eligible].n.sum() / max(len(S), 1)),
                            n_mbj_dedup_within_subset=len(Sm), K_eligible_mbj=int((nm >= min_n).sum()),
                            n_mbj_in_eligible_groups=int(G[G.eligible].n_mbj.sum()),
                            median_largest_chemsys_share_eligible=float(G[G.eligible].largest_chemsys_share.median())
                            if K else np.nan))
        for t in a["size_grid"]:
            th_rows.append(dict(release=release, scope=scope, min_group_size=int(t), K=int((G.n >= t).sum()),
                                share_of_entries=float(G[G.n >= t].n.sum() / max(len(S), 1))))
        cs = S.groupby("crys").size().reindex(CRYSTAL_SYSTEMS).fillna(0).astype(int)
        for c, n in cs.items():
            cs_rows.append(dict(release=release, scope=scope, crystal_system=c, n=int(n)))
    return (pd.concat(sp_rows, ignore_index=True), pd.DataFrame(th_rows), pd.DataFrame(cs_rows),
            pd.DataFrame(sc_rows), eligible)


def decide(SC: pd.DataFrame, a: dict, release: str) -> dict:
    s = SC[SC.release == release].set_index("scope")
    need, mincal = int(a["min_groups_for_oxides"]), int(a["min_calibration_groups"])
    k_ox, k_el, k_all = (int(s.loc[x, "K_eligible"]) for x in ("oxides_oxidation_rule", "oxides_element_rule",
                                                                 "all_chemistries"))
    primary_ox = (k_ox >= need and k_el >= need) if a.get("decision_requires_both_oxide_rules", True) else k_ox >= need
    primary = "oxides_oxidation_rule" if primary_ox else "all_chemistries"
    secondary = "all_chemistries" if primary_ox else "oxides_oxidation_rule"

    def h3(K):
        return ("valid group-level split (>= {c} calibration groups + {t} training groups)".format(c=mincal, t=K - mincal)
                if K >= 2 * mincal else "group-level jackknife+ only (Lee et al. 2023; fewer than "
                f"{2 * mincal} eligible groups)")
    borderline = (k_ox >= need) != (k_el >= need)
    return dict(release=release, rule=f"oxides primary only if K >= {need} under "
                + ("both oxide rules" if a.get("decision_requires_both_oxide_rules", True) else "the primary oxide rule"),
                K=dict(oxides_oxidation_rule=k_ox, oxides_element_rule=k_el, all_chemistries=k_all),
                oxide_definition_flips_the_decision=bool(borderline),
                primary_scope=primary, secondary_scope=secondary,
                K_primary=int(s.loc[primary, "K_eligible"]), K_secondary=int(s.loc[secondary, "K_eligible"]),
                h3_design_primary=h3(int(s.loc[primary, "K_eligible"])),
                h3_design_secondary=h3(int(s.loc[secondary, "K_eligible"])),
                tbmbj_primary=dict(n=int(s.loc[primary, "n_mbj_dedup_within_subset"]),
                                   K_eligible=int(s.loc[primary, "K_eligible_mbj"])))


# ----------------------------------------------------------------------------- main
def run_audit6(cfg: dict):
    t0 = time.time()
    import pymatgen.core
    paths = get_paths(cfg)
    outdir = os.path.join(paths.out, "audit")
    os.makedirs(outdir, exist_ok=True)
    J, a = cfg["jarvis"], cfg["audit"]
    inputs_dir = os.path.join(paths.out, J.get("inputs_subdir", "inputs"))
    releases = [J["primary"]] + list(J.get("compare", []))
    info, SP, TH, CS, SC, ST, elig = [], [], [], [], [], [], {}
    for rel in releases:
        inf = fetch_release(rel, inputs_dir)
        print(f"[audit6] {rel['name']}: {inf['bytes'] / 1e6:.1f} MB, sha256 {inf['sha256'][:12]}", flush=True)
        entries = load_release(inf)
        inf["n_entries"] = len(entries)
        info.append(inf)
        T = entry_table(entries)
        del entries
        if T.jid.duplicated().any():
            raise AssertionError(f"{rel['name']}: duplicated JARVIS ids")
        T = classify_oxides(T, a["other_anion_formers"])
        st = T[T.has_O].oxide_status.where(~T[T.has_O].oxide_status.str.startswith("other_anion"), "other_anion")
        ST.append(st.value_counts().rename_axis("oxide_status").reset_index(name="n").assign(release=rel["name"]))
        keep = ["jid", "reduced_formula", "chemsys", "n_elements", "nat", "spg_number", "spg_symbol", "crys",
                "dimensionality", "has_O", "oxide_status", "oxide_oxidation_rule", "oxide_element_rule",
                "formation_energy_peratom", "optb88vdw_total_energy", "ehull"]
        slim = T[keep].assign(has_optb88vdw_gap=T.optb88vdw_bandgap.notna(), has_mbj_gap=T.mbj_bandgap.notna())
        atomic_write_bytes(os.path.join(outdir, f"jarvis_{rel['name']}_entries.csv.gz"), _gz(slim))
        sp, th, cs, sc, el = audit_release(T, a, rel["name"])
        SP.append(sp), TH.append(th), CS.append(cs), SC.append(sc)
        elig[rel["name"]] = el
        from pymatgen.core import Composition
        uf = {f: Composition(f).reduced_formula for f in T.formula_field.dropna().unique()}
        mism = int((T.formula_field.notna() & (T.formula_field.map(uf) != T.reduced_formula)).sum())
        inf["n_formula_field_disagrees_with_atoms"] = mism
        print(f"[audit6] {rel['name']}: {len(T)} entries; 'formula' field disagrees with the atoms in {mism}; "
              f"O-containing compounds {int(T.has_O.sum())}", flush=True)
    SP, TH, CS, SC, ST = (pd.concat(x, ignore_index=True) for x in (SP, TH, CS, SC, ST))
    prim = J["primary"]["name"]
    dec = decide(SC, a, prim)
    dec_cmp = {r["name"]: decide(SC, a, r["name"]) for r in J.get("compare", [])}
    E = elig[prim][dec["primary_scope"]][["spg_number", "spg_symbol", "crys", "n", "n_mbj", "n_chemsys",
                                         "largest_chemsys_share", "share_formula_in_other_spg"]]
    atomic_write_csv(os.path.join(outdir, "audit_spacegroups.csv"), SP)
    atomic_write_csv(os.path.join(outdir, "audit_thresholds.csv"), TH)
    atomic_write_csv(os.path.join(outdir, "audit_crystal_systems.csv"), CS)
    atomic_write_csv(os.path.join(outdir, "audit_scopes.csv"), SC)
    atomic_write_csv(os.path.join(outdir, "audit_oxide_status.csv"), ST)
    atomic_write_csv(os.path.join(outdir, "audit_eligible_groups_primary.csv"), E)
    try:
        import jarvis.db.figshare as jf
        jt = {k: jf.get_db_info()[k][0] for k in ("dft_3d", "dft_3d_2022")}
    except Exception:
        jt = "jarvis-tools not installed (not needed: the figshare files are read directly)"
    meta = dict(created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), config_hash=config_hash(cfg),
                code_hash=code_hash(), pymatgen_version=getattr(pymatgen.core, "__version__", "unknown"),
                jarvis_tools_urls=jt, releases=info, audit_settings=a, decision=dec, decision_compare_releases=dec_cmp,
                note="Counts only. No target value is summarised; formation energy is used only to keep the "
                     "lowest-energy duplicate. The decision rule was fixed in the config before this run.",
                minutes=(time.time() - t0) / 60)
    atomic_write_json(os.path.join(outdir, "AUDIT_DECISION.json"), meta)

    s = SC[SC.release == prim].set_index("scope")
    lines = [f"# Phase 6 Task 6.1: JARVIS-DFT data audit ({prim}, {info[0]['n_entries']} entries)", "",
             f"Release: {info[0]['file']} ({info[0]['url']}), SHA-256 {info[0]['sha256']}.", "",
             "## Scopes after de-duplication by (reduced formula, space group)", "",
             SC.round(3).to_markdown(index=False), "",
             f"## Eligible space groups (n >= {a['min_group_size']}) by threshold", "",
             TH.pivot_table(index=["release", "min_group_size"], columns="scope", values="K").reset_index()
             .to_markdown(index=False), "",
             "## Oxide classification (O-containing compounds)", "",
             ST.pivot_table(index="oxide_status", columns="release", values="n").fillna(0).astype(int)
             .reset_index().to_markdown(index=False), "",
             "## Decision (rule fixed in the config before the run)", "", "```", json.dumps(dec, indent=2), "```", "",
             f"Other releases: {json.dumps({k: dict(K=v['K'], primary_scope=v['primary_scope']) for k, v in dec_cmp.items()})}",
             "", f"## Eligible space groups of the primary scope ({dec['primary_scope']})", "",
             E.round(3).to_markdown(index=False), "",
             "## Crystal systems (secondary analysis; only 7 groups)", "",
             CS[CS.release == prim].pivot_table(index="crystal_system", columns="scope", values="n")
             .reindex(CRYSTAL_SYSTEMS).reset_index().to_markdown(index=False), ""]
    txt = "\n".join(lines)
    atomic_write_text(os.path.join(outdir, "AUDIT_SUMMARY.md"), txt)
    print(txt)
    print(f"\n[audit6] written to {outdir} ({(time.time() - t0) / 60:.1f} min)")
    return meta
