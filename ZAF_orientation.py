"""Cartesian orientation relationships between two ZAF crystal phases.

Two *specific, signed* direction pairs define a full three-dimensional
orientation relationship.  One pair is also useful, but deliberately weaker:
it constrains that parallel zone-axis pair while leaving rotation about the
common direction undetermined.  A family such as <111> contains several
variants, so these helpers always operate on individual signed directions.
This module does not change the single-phase indexing algorithm.
"""

from __future__ import annotations

import itertools
import math
from typing import TYPE_CHECKING, Sequence

import numpy as np

if TYPE_CHECKING:
    from ZAF_crystals import CrystalModel


def _unit_direction(
    direction: Sequence[int],
    structure: str,
    model: CrystalModel | None,
) -> np.ndarray:
    # Import lazily so ZAF can import this module without a circular import.
    from ZAF import direction_cartesian, four_to_three_direction

    indices = tuple(int(value) for value in direction)
    if len(indices) == 4:
        if not (structure.upper() == "HCP" or (model is not None and model.is_hexagonal)):
            raise ValueError("Four-index directions require a hexagonal crystal.")
        indices = four_to_three_direction(indices)
    elif len(indices) != 3:
        raise ValueError("A direction needs three indices, or four for a hexagonal crystal.")
    if not any(indices):
        raise ValueError("A zero vector is not a crystal direction.")
    vector = direction_cartesian(indices, structure, crystal_model=model)
    length = float(np.linalg.norm(vector))
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("The crystal direction has invalid Cartesian length.")
    return np.asarray(vector, dtype=float) / length


def _three_index_direction(
    direction: Sequence[int],
    structure: str,
    model: CrystalModel | None,
) -> tuple[int, int, int]:
    """Validate a direction and return ZAF's primitive three-index form."""
    from ZAF import four_to_three_direction, reduce_miller

    indices = tuple(int(value) for value in direction)
    if len(indices) == 4:
        if not (structure.upper() == "HCP" or (model is not None and model.is_hexagonal)):
            raise ValueError("Four-index directions require a hexagonal crystal.")
        indices = four_to_three_direction(indices)
    elif len(indices) != 3:
        raise ValueError("A direction needs three indices, or four for a hexagonal crystal.")
    if not any(indices):
        raise ValueError("A zero vector is not a crystal direction.")
    return reduce_miller(indices)


def _axis_rotation(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rodrigues rotation about an already-normalized Cartesian axis."""
    x, y, z = np.asarray(axis, dtype=float)
    skew = np.array(((0.0, -z, y), (z, 0.0, -x), (-y, x, 0.0)), dtype=float)
    sine = math.sin(angle_rad)
    cosine = math.cos(angle_rad)
    return np.eye(3) + sine * skew + (1.0 - cosine) * (skew @ skew)


def rotation_from_parallel_pair(
    a_direction: Sequence[int],
    b_direction: Sequence[int],
    a_structure: str,
    b_structure: str,
    a_model: CrystalModel | None = None,
    b_model: CrystalModel | None = None,
    twist_deg: float = 0.0,
) -> np.ndarray:
    """Return one rotation mapping a signed phase-A direction to phase B.

    A single parallel pair has one free degree of freedom. ``twist_deg`` picks
    a rotation about the mapped phase-B direction; zero selects a stable,
    deterministic minimum-angle representative.  Callers must not interpret
    that arbitrary representative as a predicted in-plane orientation.
    """
    if not math.isfinite(twist_deg):
        raise ValueError("The single-pair twist must be finite.")
    source = _unit_direction(a_direction, a_structure, a_model)
    target = _unit_direction(b_direction, b_structure, b_model)
    cosine = float(np.clip(np.dot(source, target), -1.0, 1.0))
    cross = np.cross(source, target)
    sine = float(np.linalg.norm(cross))
    if sine > 1e-12:
        base = _axis_rotation(cross / sine, math.atan2(sine, cosine))
    elif cosine > 0.0:
        base = np.eye(3)
    else:
        # Antiparallel vectors admit infinitely many 180-degree axes. Pick the
        # Cartesian basis direction least parallel to the source, then form a
        # deterministic perpendicular axis.
        basis = np.eye(3)[int(np.argmin(np.abs(source)))]
        axis = np.cross(source, basis)
        axis /= np.linalg.norm(axis)
        base = _axis_rotation(axis, math.pi)
    return _axis_rotation(target, math.radians(twist_deg)) @ base


def target_zone_from_parallel_pair(
    observed_source_zone: Sequence[int],
    specified_source_direction: Sequence[int],
    specified_target_direction: Sequence[int],
    source_structure: str,
    target_structure: str,
    source_model: CrystalModel | None = None,
    target_model: CrystalModel | None = None,
    tolerance_deg: float = 2.0,
) -> tuple[tuple[int, int, int] | None, float]:
    """Apply a one-pair OR when the observed source is that specified axis.

    Returns the corresponding signed target zone and the source-axis mismatch.
    If the indexed source is not parallel to the entered source direction, the
    target zone is ``None``: rotation about the single common direction is free,
    so another target axis cannot be predicted uniquely.  This explicit result
    prevents callers from launching an ill-posed or unbounded twist search.
    """
    if not math.isfinite(tolerance_deg) or not 0.0 <= tolerance_deg <= 30.0:
        raise ValueError("The parallel-axis tolerance must be between 0 and 30 degrees.")
    observed = _unit_direction(observed_source_zone, source_structure, source_model)
    specified = _unit_direction(
        specified_source_direction, source_structure, source_model
    )
    # Validate the target metric as well as its indices.
    _unit_direction(specified_target_direction, target_structure, target_model)
    target_indices = _three_index_direction(
        specified_target_direction, target_structure, target_model
    )
    signed_cosine = float(np.clip(np.dot(observed, specified), -1.0, 1.0))
    absolute_cosine = abs(signed_cosine)
    line_error = (
        0.0
        if 1.0 - absolute_cosine <= 1e-12
        else math.degrees(math.acos(absolute_cosine))
    )
    if line_error > tolerance_deg:
        return None, line_error
    sign = 1 if signed_cosine >= 0.0 else -1
    return tuple(sign * value for value in target_indices), line_error


def _orientation_frame(primary: np.ndarray, secondary: np.ndarray) -> np.ndarray:
    across = secondary - float(np.dot(primary, secondary)) * primary
    across_length = float(np.linalg.norm(across))
    if across_length < 1e-4:
        raise ValueError("The two directions in each phase must be non-collinear.")
    across /= across_length
    return np.column_stack((primary, across, np.cross(primary, across)))


def rotation_from_parallel_pairs(
    a_primary: Sequence[int],
    a_secondary: Sequence[int],
    b_primary: Sequence[int],
    b_secondary: Sequence[int],
    a_structure: str,
    b_structure: str,
    a_model: CrystalModel | None = None,
    b_model: CrystalModel | None = None,
) -> np.ndarray:
    """Return a proper rotation mapping phase-A directions to phase-B directions.

    Inputs are individual signed directions, not symmetry-family labels.  The
    angle between the two phase-A directions must agree with the corresponding
    phase-B angle to within two degrees.  This tolerates a small measured or
    reported OR deviation but rejects geometrically impossible pairings.
    Four-index Miller--Bravais directions are accepted for hexagonal phases.
    """
    a1 = _unit_direction(a_primary, a_structure, a_model)
    a2 = _unit_direction(a_secondary, a_structure, a_model)
    b1 = _unit_direction(b_primary, b_structure, b_model)
    b2 = _unit_direction(b_secondary, b_structure, b_model)
    a_frame = _orientation_frame(a1, a2)
    b_frame = _orientation_frame(b1, b2)

    a_angle = math.degrees(math.acos(float(np.clip(np.dot(a1, a2), -1.0, 1.0))))
    b_angle = math.degrees(math.acos(float(np.clip(np.dot(b1, b2), -1.0, 1.0))))
    if abs(a_angle - b_angle) > 2.0:
        raise ValueError(
            "The two parallel-direction pairs have inconsistent angles "
            f"({a_angle:.2f}° in phase A versus {b_angle:.2f}° in phase B). "
            "Choose specific symmetry variants with compatible angles."
        )
    return b_frame @ a_frame.T


def predict_target_zone(
    rotation: np.ndarray,
    source_zone: Sequence[int],
    source_structure: str,
    target_structure: str,
    source_model: CrystalModel | None = None,
    target_model: CrystalModel | None = None,
    max_index: int = 8,
) -> tuple[tuple[int, int, int], float]:
    """Approximate the rotated source zone by a bounded target Miller direction.

    Returns the signed, primitive three-index target direction and its residual
    Cartesian angular error in degrees.  A nonzero error is meaningful: because
    the two cells can have different metrics, an exact low-index target zone
    need not exist.  Callers should expose or threshold this residual instead
    of silently asserting that the prediction is exact.
    """
    from ZAF import direct_lattice_matrix

    matrix = np.asarray(rotation, dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("The orientation relationship must be a finite 3×3 matrix.")
    if not np.allclose(matrix.T @ matrix, np.eye(3), atol=1e-5) or not math.isclose(
        float(np.linalg.det(matrix)), 1.0, abs_tol=1e-5
    ):
        raise ValueError("The orientation relationship must be a proper rotation.")
    if not isinstance(max_index, int) or max_index < 1 or max_index > 32:
        raise ValueError("max_index must be an integer from 1 to 32.")

    target = matrix @ _unit_direction(source_zone, source_structure, source_model)
    direct = direct_lattice_matrix(target_structure, crystal_model=target_model)
    best_zone: tuple[int, int, int] | None = None
    best_dot = -2.0
    for zone in itertools.product(range(-max_index, max_index + 1), repeat=3):
        if zone == (0, 0, 0) or math.gcd(*zone) != 1:
            continue
        vector = direct @ np.asarray(zone, dtype=float)
        similarity = float(np.dot(vector, target) / np.linalg.norm(vector))
        if similarity > best_dot + 1e-12 or (
            abs(similarity - best_dot) <= 1e-12
            and best_zone is not None
            and (sum(map(abs, zone)), max(map(abs, zone)), zone)
            < (sum(map(abs, best_zone)), max(map(abs, best_zone)), best_zone)
        ):
            best_zone, best_dot = zone, similarity
    if best_zone is None:
        raise ValueError("No target direction could be generated.")
    angular_error = math.degrees(math.acos(float(np.clip(best_dot, -1.0, 1.0))))
    return best_zone, angular_error
