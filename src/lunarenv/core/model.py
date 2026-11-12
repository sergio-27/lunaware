"""What every environment model must declare, and the registry of models.

A model is a deterministic computation with a name and a version. Besides computing, it
must be able to *describe itself*: what it assumes, what it simplifies, where it has
been validated, which parameters it takes and how uncertain they are, and how large
its numerical error is. ``describe_models()`` returns those descriptions; they are
written into the metadata of every result.

Adding a model (radiation, seismic, structural loads) means: a class with these
attributes, a parameter block in the scenario schema, and a registration call. The
runner executes models in dependency order and hands each one the outputs of those it
requires.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from pydantic import BaseModel, ConfigDict, Field


class ParameterInfo(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    default: float | str | int | None
    unit: str = ""
    description: str = ""
    uncertainty: str = Field(
        "", description="Interval or statement used for uncertainty propagation."
    )


class ModelDescription(BaseModel):
    """Self-description of a model, as published with results."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    version: str
    summary: str
    requires: tuple[str, ...] = ()
    outputs: tuple[str, ...] = ()
    assumptions: tuple[str, ...] = ()
    simplifications: tuple[str, ...] = ()
    validity_range: tuple[str, ...] = ()
    parameters: tuple[ParameterInfo, ...] = ()
    numerical_error: dict[str, float | str] = Field(default_factory=dict)
    uncertainty_sources: tuple[str, ...] = ()
    validation_cases: tuple[str, ...] = ()
    references: tuple[str, ...] = ()
    documentation: str = ""


@runtime_checkable
class EnvironmentModel(Protocol):
    """Interface of a model as seen by the registry and the runner."""

    name: str
    version: str
    requires: tuple[str, ...]

    def describe(self) -> ModelDescription: ...


_REGISTRY: dict[str, EnvironmentModel] = {}


def register(model: EnvironmentModel) -> None:
    """Add a model to the registry. Names are unique."""
    if model.name in _REGISTRY and _REGISTRY[model.name] is not model:
        raise ValueError(f"a model named {model.name!r} is already registered")
    _REGISTRY[model.name] = model


def registered_models() -> dict[str, EnvironmentModel]:
    """Registered models in an order that satisfies their dependencies."""
    import lunarenv.runner  # noqa: F401  (registers the built-in models on first use)

    ordered: dict[str, EnvironmentModel] = {}
    remaining = dict(_REGISTRY)
    while remaining:
        ready = [n for n, m in remaining.items() if all(r in ordered for r in m.requires)]
        if not ready:
            raise ValueError(f"unsatisfied or circular model dependencies: {sorted(remaining)}")
        for name in sorted(ready):
            ordered[name] = remaining.pop(name)
    return ordered
