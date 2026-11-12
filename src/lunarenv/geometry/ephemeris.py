"""Positions of the Sun and the Earth in the Moon body-fixed frame, from SPICE.

Frame choice (PLAN.md decision D6). The LOLA elevation models are expressed in the
mean-Earth/polar-axis frame of DE421. The current ephemeris is DE440, whose lunar
orientation defines a slightly different mean-Earth frame. NAIF provides
``MOON_ME_DE440_ME421``: DE440 orientation, rotated by a constant offset so that it
matches the DE421 mean-Earth frame. Using it keeps terrain and ephemeris consistent.
The residual against a pure DE421 computation is measured in the test-suite.

Aberration model. Positions are apparent (``LT+S``): light-time corrected and with
stellar aberration, as seen by the observer. Two entry points exist:

* :meth:`Ephemeris.position_from_center` - observer at the Moon's centre. One SPICE
  call per epoch serves every site: a surface point is handled by subtracting its
  body-fixed position. The stellar-aberration rotation (v/c ~ 1e-4 rad) is then applied
  to the centre-to-body vector instead of the site-to-body vector, which leaves an
  error of (lunar radius / distance) * (v/c): at most 5.2e-7 rad for the Earth and
  1.3e-9 rad for the Sun. Smaller terms (6 ms light-time across one lunar radius, site
  rotation below 5 m/s) contribute under 2e-8 rad. The measured error is asserted
  against these bounds in the test-suite.
* :meth:`Ephemeris.position_from_site` - observer at the surface point (``spkcpo``).
  Exact within the SPICE model; used as the reference in tests.

CSPICE keeps global state and is not thread-safe. Kernels are loaded through
:class:`SpiceSession`, which allows one kernel set at a time per process.
"""

from __future__ import annotations

import threading
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from types import TracebackType
from typing import Literal

import numpy as np
import spiceypy
from numpy.typing import ArrayLike, NDArray

from lunarenv.data.store import DataRecord, DataStore
from lunarenv.geometry.projection import FloatArray

Body = Literal["SUN", "EARTH"]
_KM = 1000.0


@dataclass(frozen=True, slots=True)
class KernelSet:
    """A consistent set of SPICE kernels and the lunar body-fixed frame they define."""

    name: str
    kernel_ids: tuple[str, ...]
    body_fixed_frame: str


#: Default: DE440 ephemeris and orientation, in the frame aligned to DE421 mean-Earth.
DE440 = KernelSet(
    name="DE440 / MOON_ME_DE440_ME421",
    kernel_ids=(
        "naif.lsk",
        "naif.pck.text",
        "naif.spk.de440s",
        "naif.pck.moon_pa_de440",
        "naif.fk.moon_de440",
    ),
    body_fixed_frame="MOON_ME_DE440_ME421",
)

#: The ephemeris and frame of the LOLA products themselves. Reference for decision D6.
DE421 = KernelSet(
    name="DE421 / MOON_ME_DE421",
    kernel_ids=(
        "naif.lsk",
        "naif.pck.text",
        "naif.spk.de421",
        "naif.pck.moon_pa_de421",
        "naif.fk.moon_de421",
    ),
    body_fixed_frame="MOON_ME_DE421",
)

_session_lock = threading.Lock()
_active_session: SpiceSession | None = None


class SpiceSession:
    """Context manager that owns the process-wide SPICE kernel pool.

    Entering loads exactly the kernels of ``kernel_set`` from the data cache; leaving
    unloads them. The pool must be empty on entry: kernels loaded by other code could
    silently change frames or constants, and clearing them would break that code, so
    the session refuses to start instead. Sessions cannot be nested or run concurrently,
    because the kernel pool is global.
    """

    def __init__(self, store: DataStore, kernel_set: KernelSet = DE440) -> None:
        self.store = store
        self.kernel_set = kernel_set
        self._loaded: list[str] = []

    def __enter__(self) -> SpiceSession:
        global _active_session
        with _session_lock:
            if _active_session is not None:
                raise RuntimeError("a SPICE session is already active in this process")
            if spiceypy.ktotal("ALL") != 0:
                raise RuntimeError(
                    "SPICE kernels are already loaded in this process; lunarenv needs an "
                    "empty kernel pool to guarantee which frames and constants are used"
                )
            paths = [str(self.store.path(kernel_id)) for kernel_id in self.kernel_set.kernel_ids]
            try:
                for path in paths:
                    spiceypy.furnsh(path)
                    self._loaded.append(path)
            except BaseException:
                self._unload()
                raise
            _active_session = self
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        global _active_session
        with _session_lock:
            self._unload()
            _active_session = None

    def _unload(self) -> None:
        while self._loaded:
            spiceypy.unload(self._loaded.pop())

    @property
    def active(self) -> bool:
        return bool(self._loaded)

    def provenance(self) -> list[DataRecord]:
        """Manifest records of the loaded kernels."""
        return self.store.provenance(self.kernel_set.kernel_ids)

    def toolkit_version(self) -> str:
        return str(spiceypy.tkvrsn("TOOLKIT"))


@contextmanager
def spice_session(store: DataStore, kernel_set: KernelSet = DE440) -> Iterator[Ephemeris]:
    """Open a SPICE session and yield an :class:`Ephemeris` bound to it."""
    with SpiceSession(store, kernel_set) as session:
        yield Ephemeris(session)


class Ephemeris:
    """Apparent positions of the Sun and the Earth in the lunar body-fixed frame."""

    def __init__(self, session: SpiceSession) -> None:
        if not session.active:
            raise RuntimeError("the SPICE session is not active")
        self._session = session
        self.frame = session.kernel_set.body_fixed_frame

    def _require_active(self) -> None:
        if not self._session.active:
            raise RuntimeError("the SPICE session that created this object has been closed")

    # ---------------------------------------------------------------------- time
    def utc_to_et(self, utc: ArrayLike) -> FloatArray:
        """UTC to ephemeris time (TDB seconds past J2000).

        Accepts ISO-8601 strings or ``numpy.datetime64`` values, scalar or array.
        Leap seconds come from the loaded LSK.
        """
        self._require_active()
        values = np.asarray(utc)
        if values.dtype.kind == "M":
            strings = np.datetime_as_string(values, unit="us")
        elif values.dtype.kind in ("U", "S", "O"):
            strings = values.astype(str)
        else:
            raise TypeError("utc must be ISO-8601 strings or numpy datetime64")
        flat = [float(spiceypy.utc2et(s)) for s in np.atleast_1d(strings).ravel()]
        return np.asarray(flat, dtype=np.float64).reshape(values.shape)

    def et_to_utc(self, et: ArrayLike) -> NDArray[np.str_]:
        """Ephemeris time to ISO-8601 UTC strings with millisecond resolution."""
        self._require_active()
        values = np.asarray(et, dtype=np.float64)
        flat = [str(spiceypy.et2utc(float(t), "ISOC", 3)) for t in np.atleast_1d(values).ravel()]
        return np.asarray(flat).reshape(values.shape)

    # ----------------------------------------------------------------- positions
    def position_from_center(self, body: Body, et: ArrayLike) -> FloatArray:
        """Apparent position of ``body`` relative to the Moon's centre, metres.

        Shape ``(..., 3)`` in the body-fixed frame evaluated at the observation epoch.
        """
        self._require_active()
        epochs = np.asarray(et, dtype=np.float64)
        out = np.empty((epochs.size, 3), dtype=np.float64)
        for i, t in enumerate(epochs.ravel()):
            position_km, _ = spiceypy.spkpos(body, float(t), self.frame, "LT+S", "MOON")
            out[i] = position_km
        return (out * _KM).reshape(*epochs.shape, 3)

    def position_from_site(
        self, body: Body, et: ArrayLike, site_position_m: Sequence[float] | FloatArray
    ) -> FloatArray:
        """Apparent position of ``body`` relative to a surface point, metres.

        ``site_position_m`` is the body-fixed Cartesian position of the observer. The
        observer's own motion (lunar rotation and orbit) enters the aberration
        correction. Shape ``(..., 3)``.
        """
        self._require_active()
        site_km = np.asarray(site_position_m, dtype=np.float64) / _KM
        if site_km.shape != (3,):
            raise ValueError("site position must be a 3-vector")
        epochs = np.asarray(et, dtype=np.float64)
        out = np.empty((epochs.size, 3), dtype=np.float64)
        for i, t in enumerate(epochs.ravel()):
            state_km, _ = spiceypy.spkcpo(
                body, float(t), self.frame, "OBSERVER", "LT+S", site_km, "MOON", self.frame
            )
            out[i] = state_km[:3]
        return (out * _KM).reshape(*epochs.shape, 3)

    # --------------------------------------------------------------------- sizes
    def body_radius_m(self, body: Body | Literal["MOON"]) -> float:
        """Mean of the tri-axial radii in the loaded PCK, metres.

        For the Sun the PCK value is the photospheric radius used to define the visible
        disc (695 700 km in pck00011, the IAU 2015 nominal value).
        """
        self._require_active()
        _, radii_km = spiceypy.bodvrd(body, "RADII", 3)
        return float(np.mean(radii_km)) * _KM

    def frame_rotation(self, from_frame: str, et: float) -> FloatArray:
        """Rotation matrix taking vectors from ``from_frame`` to the body-fixed frame."""
        self._require_active()
        return np.asarray(spiceypy.pxform(from_frame, self.frame, float(et)), dtype=np.float64)
