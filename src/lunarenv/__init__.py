"""lunarenv: open, verified lunar surface environment engine.

Public entry points (imported lazily, so that ``import lunarenv`` stays cheap):
``run_scenario``, ``describe_models``, ``get_validation_status``.
"""

from importlib.metadata import PackageNotFoundError, version
from typing import Any

try:
    __version__ = version("lunarenv")
except PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "0.0.0+unknown"

_API = ("run_scenario", "describe_models", "get_validation_status", "load_scenario", "load_run")
__all__ = ["__version__", *_API]


def __getattr__(name: str) -> Any:
    if name in _API:
        from lunarenv import api

        return getattr(api, name)
    raise AttributeError(f"module 'lunarenv' has no attribute {name!r}")
