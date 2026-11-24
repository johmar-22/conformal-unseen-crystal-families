"""Configuration loading, output paths, hashing and small I/O helpers."""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass

import yaml

PKG_DIR = os.path.dirname(os.path.abspath(__file__))
BUNDLE_DIR = os.path.dirname(PKG_DIR)
DEFAULT_CONFIG = os.path.join(BUNDLE_DIR, "config", "phase05.yaml")


def load_config(path: str | None = None) -> dict:
    path = path or DEFAULT_CONFIG
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg["_config_path"] = os.path.abspath(path)
    return cfg


def config_hash(cfg: dict) -> str:
    """Hash of the scientific settings (paths and compute settings excluded)."""
    keep = {k: v for k, v in cfg.items() if k not in ("paths", "compute") and not k.startswith("_")}
    return hashlib.sha256(json.dumps(keep, sort_keys=True).encode()).hexdigest()


def code_hash() -> str:
    """Hash of every .py file in the package, so results can be tied to the exact code."""
    h = hashlib.sha256()
    for name in sorted(os.listdir(PKG_DIR)):
        if name.endswith(".py"):
            with open(os.path.join(PKG_DIR, name), "rb") as f:
                h.update(name.encode())
                h.update(f.read())
    return h.hexdigest()


def sha256_file(path: str, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


@dataclass
class Paths:
    drive_root: str
    out: str
    cache: str
    phase2: str            # Phase 2 output folder (read only): LOSO units, de-duplication record

    @property
    def inputs(self): return self._mk("inputs")
    @property
    def phase0(self): return self._mk("phase0")
    @property
    def scores(self): return self._mk("scores")
    @property
    def report(self): return self._mk("report")

    @property
    def phase2_loso(self): return os.path.join(self.phase2, "loso")

    def _mk(self, sub: str) -> str:
        d = os.path.join(self.out, sub)
        os.makedirs(d, exist_ok=True)
        return d

    def drive_file(self, name: str) -> str:
        return os.path.join(self.drive_root, name)


def get_paths(cfg: dict) -> Paths:
    p = cfg["paths"]
    root = p["drive_root"]
    if not os.path.isdir(root):
        raise FileNotFoundError(
            f"Drive folder not found: {root}\nMount Drive first and check paths.drive_root in the config.")
    out = os.path.join(root, p["out_dir"])
    os.makedirs(out, exist_ok=True)
    os.makedirs(p["local_cache"], exist_ok=True)
    return Paths(drive_root=root, out=out, cache=p["local_cache"], phase2=os.path.join(root, p["phase2_dir"]))


def atomic_write_bytes(path: str, data: bytes) -> None:
    """Write via a temp file and rename, so a Colab disconnect never leaves a half-written file."""
    d = os.path.dirname(path)
    os.makedirs(d, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=d, suffix=".tmp")
    with os.fdopen(fd, "wb") as f:
        f.write(data)
    os.replace(tmp, path)


def atomic_write_json(path: str, obj) -> None:
    atomic_write_bytes(path, json.dumps(obj, indent=2, default=_json_default).encode())


def atomic_write_text(path: str, text: str) -> None:
    atomic_write_bytes(path, text.encode())


def atomic_write_csv(path: str, df) -> None:
    atomic_write_bytes(path, df.to_csv(index=False).encode())


def atomic_save_npz(path: str, **arrays) -> None:
    import io
    import numpy as np
    buf = io.BytesIO()
    np.savez_compressed(buf, **arrays)
    atomic_write_bytes(path, buf.getvalue())


def _json_default(o):
    import numpy as np
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)
