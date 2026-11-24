"""LightGBM fitting (parallel, memory-mapped features) and Optuna hyperparameter search.

Fixes relative to the v5 pipeline
---------------------------------
* v5 passed `subsample` without `subsample_freq`. LightGBM only bags when bagging_freq > 0,
  so the reported subsample ratios had no effect. Here bagging_freq = 1, so the tuned
  bagging fraction is actually used.
* v5 tuned with early stopping but refit with a fixed 700 rounds. Here the refit uses the
  mean best iteration of the selected trial's inner folds.
* Hyperparameters are tuned on training data only (never on calibration or test data).
  v2: the inner folds are passed in explicitly; in Phase 2 they are leave-one-lattice-out folds
  over the five training lattices (nested LOSO).
"""
from __future__ import annotations

import json
import os
import warnings

import numpy as np

from .config import atomic_write_json

# loky restarts idle workers as a memory-leak guard and warns each time; results are unaffected
warnings.filterwarnings("ignore", message="A worker stopped while some jobs", category=UserWarning)

SEARCH_SPACE = {  # identical ranges to v5 (manuscript Table 1)
    "learning_rate": ("float", 0.01, 0.15, True),
    "num_leaves": ("int", 31, 127, False),
    "max_depth": ("int", 5, 9, False),
    "min_child_samples": ("int", 20, 60, False),
    "subsample": ("float", 0.65, 0.95, False),
    "colsample_bytree": ("float", 0.65, 0.95, False),
    "reg_alpha": ("float", 1e-4, 1.0, True),
    "reg_lambda": ("float", 1e-4, 1.0, True),
}
_NATIVE = {"min_child_samples": "min_data_in_leaf", "subsample": "bagging_fraction",
           "colsample_bytree": "feature_fraction", "reg_alpha": "lambda_l1", "reg_lambda": "lambda_l2"}


def n_cpu() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except Exception:
        return os.cpu_count() or 1


def plan_parallel(cfg: dict) -> tuple[int, int]:
    """(concurrent fits, threads per fit). LightGBM on ~1.5k rows scales poorly past ~4 threads,
    so on a 12-core Colab A100 host three 4-thread fits beat one 12-thread fit."""
    c = n_cpu()
    w = cfg["compute"].get("parallel_fits", "auto")
    w = max(1, c // 4) if w in (None, "auto") else max(1, int(w))
    return w, max(1, c // w)


def lgb_params(hp: dict, cfg: dict, seed: int, threads: int, objective: str | None = None,
               quantile: float | None = None) -> dict:
    p = {"objective": objective or cfg["lgbm"]["objective"], "verbosity": -1, "seed": int(seed),
         "num_threads": int(threads), "deterministic": bool(cfg["lgbm"]["deterministic"]),
         "force_col_wise": True, "bagging_freq": 1, "bagging_seed": int(seed),
         "feature_fraction_seed": int(seed)}
    for k, v in hp.items():
        if k == "n_rounds" or k.startswith("_"):
            continue
        p[_NATIVE.get(k, k)] = v
    if quantile is not None:
        p["objective"], p["alpha"] = "quantile", float(quantile)
    if cfg["lgbm"].get("device", "cpu") == "gpu":
        p["device_type"] = "gpu"
    return p


def _dataset(X, y, cfg_max_bin):
    import lightgbm as lgb
    return lgb.Dataset(X, label=y, params={"max_bin": int(cfg_max_bin), "verbosity": -1,
                                           "feature_pre_filter": False}, free_raw_data=True)


def fit_predict_task(task: dict) -> dict:
    """Worker: fit one booster and return predictions for the requested row sets only
    (boosters are never pickled back, which keeps memory flat)."""
    import lightgbm as lgb
    X = np.load(task["X_path"], mmap_mode="r")
    dtr = _dataset(np.ascontiguousarray(X[task["rows_train"]]), task["y_train"], task["max_bin"])
    kw = {}
    if task.get("rows_valid") is not None:
        dva = _dataset(np.ascontiguousarray(X[task["rows_valid"]]), task["y_valid"], task["max_bin"])
        kw = dict(valid_sets=[dva], callbacks=[lgb.early_stopping(task["early_stop"], verbose=False)])
        task["params"]["metric"] = "l1"
    booster = lgb.train(task["params"], dtr, num_boost_round=int(task["n_rounds"]), **kw)
    out = {"best_iteration": int(booster.best_iteration or task["n_rounds"])}
    for name, rows in task.get("predict", {}).items():
        out[name] = booster.predict(np.ascontiguousarray(X[rows]), num_threads=task["params"]["num_threads"])
    return out


def run_tasks(tasks: list[dict], n_workers: int) -> list[dict]:
    if n_workers == 1 or len(tasks) == 1:
        return [fit_predict_task(t) for t in tasks]
    from joblib import Parallel, delayed
    with warnings.catch_warnings():
        # loky restarts idle workers as a memory-leak guard and warns about it; results are unaffected
        warnings.filterwarnings("ignore", message="A worker stopped while some jobs")
        return Parallel(n_jobs=n_workers, backend="loky", max_nbytes=None)(
            delayed(fit_predict_task)(t) for t in tasks)


# ----------------------------------------------------------------------------- HPO
def _distributions():
    import optuna
    d = {}
    for k, (kind, lo, hi, log) in SEARCH_SPACE.items():
        d[k] = (optuna.distributions.FloatDistribution(lo, hi, log=log) if kind == "float"
                else optuna.distributions.IntDistribution(lo, hi))
    return d


def tune(X_path: str, rows: np.ndarray, y: np.ndarray, folds: list, h: dict, cfg: dict, seed: int,
         out_json: str, log_json: str, label: str = "") -> dict:
    """Optuna TPE over SEARCH_SPACE, MAE objective averaged over the given inner folds
    (list of (train_idx, valid_idx) positions into `rows`). Completed trials are written to
    `log_json` after each trial, so a disconnect loses at most one trial."""
    import optuna

    if os.path.exists(out_json):
        with open(out_json) as f:
            return json.load(f)
    workers, threads = plan_parallel(cfg)
    dists = _distributions()
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction="minimize", sampler=optuna.samplers.TPESampler(
        seed=seed, n_startup_trials=h["n_startup_trials"]))
    done = []
    if os.path.exists(log_json):
        with open(log_json) as f:
            done = json.load(f)
        for t in done:
            study.add_trial(optuna.trial.create_trial(params=t["params"], distributions=dists,
                                                      value=t["value"], user_attrs={"iters": t["iters"]}))

    def objective(trial):
        hp = {k: (trial.suggest_float(k, lo, hi, log=log) if kind == "float" else trial.suggest_int(k, lo, hi))
              for k, (kind, lo, hi, log) in SEARCH_SPACE.items()}
        tasks = []
        for tr, va in folds:
            tasks.append(dict(X_path=X_path, max_bin=cfg["lgbm"]["max_bin"],
                              rows_train=rows[tr], y_train=y[tr], rows_valid=rows[va], y_valid=y[va],
                              early_stop=h["early_stopping_rounds"], n_rounds=h["max_rounds"],
                              params=lgb_params(hp, cfg, seed + trial.number, threads),
                              predict={"valid": rows[va]}))
        res = run_tasks(tasks, workers)
        maes = [float(np.mean(np.abs(y[va] - r["valid"]))) for (_, va), r in zip(folds, res)]
        trial.set_user_attr("iters", [r["best_iteration"] for r in res])
        return float(np.mean(maes))

    def save(study_, trial):
        if trial.state == optuna.trial.TrialState.COMPLETE:
            done.append({"params": trial.params, "value": trial.value, "iters": trial.user_attrs["iters"]})
            atomic_write_json(log_json, done)
            print(f"    {label} trial {len(done):>2}/{h['n_trials']}  MAE={trial.value:.5f}  "
                  f"best={study_.best_value:.5f}", flush=True)

    remaining = h["n_trials"] - len(done)
    if remaining > 0:
        study.optimize(objective, n_trials=remaining, callbacks=[save])
    best = study.best_trial
    result = dict(best.params)
    result["n_rounds"] = int(max(1, round(float(np.mean(best.user_attrs["iters"])))))
    result["_cv_mae"] = float(best.value)
    result["_inner_best_iters"] = best.user_attrs["iters"]
    atomic_write_json(out_json, result)
    return result
