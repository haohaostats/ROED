"""Load the frozen exact engine bundled with this reproduction archive."""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path


ENGINE_PATH = Path(__file__).resolve().with_name("exact_engine.py")


def engine_digest() -> str:
    return hashlib.sha256(ENGINE_PATH.read_bytes()).hexdigest()


def load_engine():
    name = "roed_paper_exact_engine"
    if name in sys.modules:
        return sys.modules[name]
    if not ENGINE_PATH.exists():
        raise FileNotFoundError(f"Bundled exact engine not found: {ENGINE_PATH}")
    spec = importlib.util.spec_from_file_location(name, ENGINE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load exact engine from {ENGINE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
