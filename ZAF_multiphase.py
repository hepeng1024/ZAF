"""Joint presentation and optional crystallographic constraints for phase fits.

Each phase retains its own calculated reciprocal lattice and indexing result.
An orientation relationship is a prior supplied as one or two *specific*,
signed parallel direction pairs, not an assertion that whole symmetry
families are parallel.  Two pairs determine a complete 3-D relationship; one
pair constrains only that corresponding zone and leaves the in-plane twist
free.  The image itself remains the evidence for whether either phase is
present; shared spots are shown for both phases rather than silently claimed
by one of them.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Sequence

import numpy as np
from PIL import Image, ImageDraw

import ZAF
from ZAF_orientation import (
    predict_target_zone,
    rotation_from_parallel_pair,
    rotation_from_parallel_pairs,
    target_zone_from_parallel_pair,
)

if TYPE_CHECKING:
    from ZAF_crystals import CrystalModel


@dataclass(frozen=True)
class PhaseSpec:
    key: str
    label: str
    crystal_structure: str
    crystal_model: CrystalModel | None = None
    known_zone: tuple[int, int, int] | None = None
    lattice_nm: float | None = None


@dataclass(frozen=True)
class OrientationSpec:
    source_key: str
    target_key: str
    # The source/target names are retained for API compatibility. In the GUI
    # these are presented neutrally as Phase 1 and Phase 2.
    source_primary: tuple[int, ...] | None = None
    source_secondary: tuple[int, ...] | None = None
    target_primary: tuple[int, ...] | None = None
    target_secondary: tuple[int, ...] | None = None


@dataclass
class PhaseResult:
    spec: PhaseSpec
    indexing: ZAF.IndexingResult
    fitted_image: Image.Image
    target_rows: list[dict[str, object]]
    independent_indexing: ZAF.IndexingResult | None = None
    or_predicted_zone: tuple[int, int, int] | None = None
    or_residual_deg: float | None = None
    or_inplane_error_deg: float | None = None


@dataclass
class MultiPhaseResult:
    phase_results: list[PhaseResult]
    composite_image: Image.Image
    warnings: list[str] = field(default_factory=list)


def _phase_lattice_nm(spec: PhaseSpec) -> float | None:
    if spec.lattice_nm is not None:
        return float(spec.lattice_nm)
    if spec.crystal_model is not None:
        return float(spec.crystal_model.a_nm)
    return None


def _phase_hcp_notation(spec: PhaseSpec) -> bool:
    return ZAF.normalize_crystal_structure(spec.crystal_structure) == "HCP"


def _line_angle_error_deg(first: float, second: float) -> float:
    """Smallest difference of two unoriented diffraction-line angles."""
    return abs((first - second + 90.0) % 180.0 - 90.0)


def _index_phase(
    gray: np.ndarray,
    image_size: tuple[int, int],
    spec: PhaseSpec,
    *,
    current_zone: tuple[int, int, int] | None,
    n_peaks: int,
    peak_percentile: float,
    scale_bar_pixels: float | None,
    scale_bar_value_inv_nm: float | None,
) -> ZAF.IndexingResult:
    return ZAF.index_diffraction_pattern(
        gray,
        image_size,
        current_zone=current_zone,
        center_xy=None,
        n_peaks=n_peaks,
        min_distance_px=None,
        spot_sigma_px=None,
        peak_percentile=peak_percentile,
        max_index=ZAF.DEFAULT_REFLECTION_MAX_INDEX,
        max_g_norm=ZAF.DEFAULT_REFLECTION_MAX_G_NORM,
        tolerance_fraction=0.18,
        scale_bar_pixels=scale_bar_pixels,
        scale_bar_value_inv_nm=scale_bar_value_inv_nm,
        expected_lattice_parameter_nm=_phase_lattice_nm(spec),
        crystal_structure=spec.crystal_structure,
        hcp_four_index=_phase_hcp_notation(spec),
        crystal_model=spec.crystal_model,
    )


def _project_or_fit(
    source_match: ZAF.PatternMatch,
    source_spec: PhaseSpec,
    target_indexing: ZAF.IndexingResult,
    target_spec: PhaseSpec,
    rotation: np.ndarray,
    image_size: tuple[int, int],
    alpha_deg: float,
    beta_deg: float,
    image_to_holder_rotation_deg: float,
    holder_order: str,
) -> tuple[ZAF.PatternMatch, float]:
    """Project the source image frame into the target lattice via the OR.

    This uses both parallel direction pairs to predict the *in-plane* target
    spot orientation, not merely its zone.  A small similarity refinement is
    accepted only if it stays close to that prior.
    """
    reference = target_indexing.best
    _zero_holder, source_image_axes, _holder_axes = ZAF.orientation_from_match(
        source_match,
        alpha_deg,
        beta_deg,
        image_to_holder_rotation_deg,
        holder_order,
    )
    target_axes = rotation @ source_image_axes
    bx = reference.basis_x
    by = reference.basis_y
    predicted_linear = np.asarray(
        [
            [float(np.dot(bx, target_axes[:, 0])), float(np.dot(by, target_axes[:, 0]))],
            [float(np.dot(bx, target_axes[:, 1])), float(np.dot(by, target_axes[:, 1]))],
        ],
        dtype=float,
    )
    u, _singular_values, vt = np.linalg.svd(predicted_linear)
    predicted_linear = u @ vt
    if float(np.linalg.det(predicted_linear)) < 0:
        u[:, -1] *= -1
        predicted_linear = u @ vt
    predicted_angle = math.degrees(
        math.atan2(predicted_linear[1, 0], predicted_linear[0, 0])
    )
    source_a = _phase_lattice_nm(source_spec)
    target_a = _phase_lattice_nm(target_spec)
    expected_scale = (
        source_match.scale * source_a / target_a
        if source_a is not None and target_a is not None
        else reference.scale
    )
    peaks_screen = ZAF.as_screen_points(target_indexing.peaks)
    fixed = ZAF.score_transform(
        reference.reflections,
        peaks_screen,
        predicted_linear,
        expected_scale,
        source_match.translation,
        image_size,
        0.18,
    )
    score, matched, visible, rms, tolerance = fixed
    best = replace(
        reference,
        linear_unit=predicted_linear,
        scale=expected_scale,
        translation=source_match.translation.copy(),
        score=score,
        matched_indices=matched,
        visible_indices=visible,
        rms_px=rms,
        tolerance_px=tolerance,
        estimated_lattice_nm=None,
        lattice_relative_error=None,
        selection_score=None,
    )
    refined = ZAF.refine_match(
        reference.reflections,
        peaks_screen,
        predicted_linear,
        expected_scale,
        source_match.translation,
        image_size,
        0.18,
    )
    unit, scale, center, score, matched, visible, rms, tolerance = refined
    refined_angle = math.degrees(math.atan2(unit[1, 0], unit[0, 0]))
    if (
        score >= best.score
        and _line_angle_error_deg(refined_angle, predicted_angle) <= 5.0
        and abs(math.log(scale / expected_scale)) <= math.log(1.08)
        and float(np.linalg.norm(center - source_match.translation)) <= 0.02 * min(image_size)
    ):
        best = replace(
            best,
            linear_unit=unit,
            scale=scale,
            translation=center,
            score=score,
            matched_indices=matched,
            visible_indices=visible,
            rms_px=rms,
            tolerance_px=tolerance,
        )
    return best, predicted_angle


def _render_composite(
    image: Image.Image, phase_results: Sequence[PhaseResult]
) -> Image.Image:
    out = image.convert("RGB").copy()
    draw = ImageDraw.Draw(out)
    colors = (
        (0, 230, 255),
        (255, 90, 190),
        (255, 210, 50),
        (110, 235, 105),
        (190, 145, 255),
        (255, 135, 80),
    )
    size = min(out.size)
    radius = max(6, int(size / 260))
    line_width = max(2, radius // 3)
    font = ZAF.load_label_font(max(22, int(size / 82)))
    label_height = max(30, int(size / 55))
    legend_width = min(out.width - 10, max(360, int(size * 0.5)))
    draw.rectangle(
        (6, 6, legend_width, 24 + label_height * (len(phase_results) + 1)),
        fill=(0, 0, 0),
    )
    draw.text((15, 12), "Matched diffraction spots by phase", font=font, fill=(255, 255, 255))
    for phase_number, phase in enumerate(phase_results):
        color = colors[phase_number % len(colors)]
        match = phase.indexing.best
        pred = ZAF.transformed_pattern_points(match, out.size, rotated=True)
        for reflection_index in match.matched_indices:
            x = float(pred[reflection_index, 0])
            y = float(-pred[reflection_index, 1])
            if not (-radius <= x <= out.width + radius and -radius <= y <= out.height + radius):
                continue
            # Concentric open rings preserve the underlying image and make a
            # spot assigned to multiple phases visibly shared, not exclusive.
            ring_radius = radius + 3 * phase_number
            draw.ellipse(
                (x - ring_radius, y - ring_radius, x + ring_radius, y + ring_radius),
                outline=color,
                width=line_width,
            )
        zone_label = ZAF.format_zone_direction(
            match.zone, match.crystal_structure, match.hcp_four_index
        )
        draw.text(
            (15, 12 + label_height * (phase_number + 1)),
            f"{phase.spec.label}: {zone_label} ({match.matched_count} matched)",
            font=font,
            fill=color,
        )
    return out


def analyze_multiphase(
    gray: np.ndarray,
    rgb: Image.Image,
    phases: Sequence[PhaseSpec],
    *,
    n_peaks: int = 80,
    peak_percentile: float = 99.0,
    scale_bar_pixels: float | None = None,
    scale_bar_value_inv_nm: float | None = None,
    orientation: OrientationSpec | None = None,
    alpha_deg: float = 0.0,
    beta_deg: float = 0.0,
    alpha_limits: tuple[float, float] = (-35.0, 35.0),
    beta_limits: tuple[float, float] = (-20.0, 20.0),
    image_to_holder_rotation_deg: float = 90.0,
    holder_order: str = "xy",
) -> MultiPhaseResult:
    """Fit selected phases to one image and report per-phase target axes.

    A supplied OR can propose a zone and re-fit that phase. With two direction
    pairs it can also constrain the in-plane variant. With one pair, only the
    entered parallel-axis correspondence is applied; twist about that axis is
    intentionally left to the image fit. It does not manufacture evidence: if
    the constrained diffraction fit is much worse than the independent fit,
    the independent result is retained and a warning explains the conflict.
    """
    if len(phases) < 2:
        raise ValueError("Multiphase analysis requires at least two phases.")
    if len({phase.key for phase in phases}) != len(phases):
        raise ValueError("Each selected phase must have a unique key.")
    if tuple(gray.shape[::-1]) != tuple(rgb.size):
        raise ValueError("The grayscale and color diffraction images have different sizes.")
    for phase in phases:
        if not phase.key.strip() or not phase.label.strip():
            raise ValueError("Every selected phase needs a key and a display label.")
        if phase.lattice_nm is not None and (not math.isfinite(phase.lattice_nm) or phase.lattice_nm <= 0):
            raise ValueError(f"{phase.label}: lattice a must be positive and finite.")

    # Validate OR completeness before any costly image indexing. In
    # particular, a deliberately blank pair 2 means a valid one-pair OR; one
    # half-filled pair 2 is an immediate input error rather than a search.
    specs = {phase.key: phase for phase in phases}
    if orientation is not None:
        if orientation.source_key not in specs or orientation.target_key not in specs:
            raise ValueError("Both orientation-relationship phases must be selected.")
        if orientation.source_key == orientation.target_key:
            raise ValueError("The orientation relationship requires two distinct phases.")
        if orientation.source_primary is None or orientation.target_primary is None:
            raise ValueError(
                "Orientation-relationship pair 1 needs one signed direction for each phase."
            )
        if (orientation.source_secondary is None) != (
            orientation.target_secondary is None
        ):
            raise ValueError(
                "Orientation-relationship pair 2 must contain directions for both "
                "phases, or be left blank for both phases."
            )

    warnings: list[str] = []
    independent: dict[str, ZAF.IndexingResult] = {}
    for phase in phases:
        try:
            independent[phase.key] = _index_phase(
                gray,
                rgb.size,
                phase,
                current_zone=phase.known_zone,
                n_peaks=n_peaks,
                peak_percentile=peak_percentile,
                scale_bar_pixels=scale_bar_pixels,
                scale_bar_value_inv_nm=scale_bar_value_inv_nm,
            )
        except (RuntimeError, ValueError) as exc:
            warnings.append(f"{phase.label}: could not fit this phase ({exc}).")
    if not independent:
        raise RuntimeError("None of the selected phases could be fitted to the diffraction image.")

    selected = dict(independent)
    predicted_zone: dict[str, tuple[tuple[int, int, int], float]] = {}
    inplane_error: dict[str, float] = {}

    def constrain_phase_from_axis(
        anchor: PhaseSpec,
        constrained_phase: PhaseSpec,
        zone: tuple[int, int, int],
    ) -> None:
        """Try one bounded image fit for an OR-predicted zone.

        Comparing this single constrained candidate with the independent fit
        makes a one-pair OR image-assisted without searching the undetermined
        continuum of twists about the parallel direction.
        """
        if (
            constrained_phase.known_zone is not None
            and ZAF.canonical_line(constrained_phase.known_zone)
            != ZAF.canonical_line(zone)
        ):
            warnings.append(
                f"{constrained_phase.label}: the entered known zone conflicts with "
                "the orientation-relationship prediction "
                f"{ZAF.format_zone_direction(zone, constrained_phase.crystal_structure, _phase_hcp_notation(constrained_phase))}; "
                "the entered zone was kept."
            )
            return
        if constrained_phase.known_zone is not None:
            return
        anchor_indexing = selected.get(anchor.key)
        if anchor_indexing is None:
            return
        anchor_match = anchor_indexing.best
        # When both lattice constants are known, the anchor fit supplies a
        # physical camera scale even when the image has no printed scale bar.
        anchor_a = _phase_lattice_nm(anchor)
        constrained_a = _phase_lattice_nm(constrained_phase)
        implied_pixels = (
            anchor_match.scale * anchor_a
            if anchor_a is not None
            and constrained_a is not None
            and scale_bar_value_inv_nm is None
            else scale_bar_pixels
        )
        implied_value = (
            1.0
            if implied_pixels is not None and scale_bar_value_inv_nm is None
            else scale_bar_value_inv_nm
        )
        try:
            constrained = _index_phase(
                gray,
                rgb.size,
                constrained_phase,
                current_zone=zone,
                n_peaks=n_peaks,
                peak_percentile=peak_percentile,
                scale_bar_pixels=implied_pixels,
                scale_bar_value_inv_nm=implied_value,
            )
        except (RuntimeError, ValueError) as exc:
            warnings.append(
                f"{constrained_phase.label}: the orientation-relationship "
                f"zone could not be fitted ({exc})."
            )
            return
        prior = selected.get(constrained_phase.key)
        if prior is None or (
            constrained.best.matched_count
            >= max(4, int(0.65 * prior.best.matched_count))
            and constrained.best.score
            >= prior.best.score - max(3.0, 0.15 * abs(prior.best.score))
        ):
            selected[constrained_phase.key] = constrained
        else:
            warnings.append(
                f"{constrained_phase.label}: diffraction spots do not support the "
                "orientation-relationship prediction "
                f"{ZAF.format_zone_direction(zone, constrained_phase.crystal_structure, _phase_hcp_notation(constrained_phase))} "
                "as well as the independent fit; check the signed direction "
                "variant and the phase overlay."
            )

    if orientation is not None:
        source = specs[orientation.source_key]
        target = specs[orientation.target_key]
        # These were validated before indexing; the assertions also narrow the
        # optional dataclass fields for static type checkers.
        assert orientation.source_primary is not None
        assert orientation.target_primary is not None
        source_has_second = orientation.source_secondary is not None
        target_has_second = orientation.target_secondary is not None
        assert source_has_second == target_has_second

        if not source_has_second:
            # Validate both directions and explicitly acknowledge the free
            # rotation about their common axis. No twist grid is searched.
            rotation_from_parallel_pair(
                orientation.source_primary,
                orientation.target_primary,
                source.crystal_structure,
                target.crystal_structure,
                source.crystal_model,
                target.crystal_model,
            )
            constrained_key: str | None = None
            source_indexing = selected.get(source.key)
            if source_indexing is not None:
                zone, error_deg = target_zone_from_parallel_pair(
                    source_indexing.best.zone,
                    orientation.source_primary,
                    orientation.target_primary,
                    source.crystal_structure,
                    target.crystal_structure,
                    source.crystal_model,
                    target.crystal_model,
                )
                if zone is not None:
                    predicted_zone[target.key] = zone, error_deg
                    constrain_phase_from_axis(source, target, zone)
                    constrained_key = target.key
            if constrained_key is None:
                # The relationship is physically symmetric. If Phase 2 is on
                # its specified axis, it can constrain Phase 1 instead.
                target_indexing = selected.get(target.key)
                if target_indexing is not None:
                    zone, error_deg = target_zone_from_parallel_pair(
                        target_indexing.best.zone,
                        orientation.target_primary,
                        orientation.source_primary,
                        target.crystal_structure,
                        source.crystal_structure,
                        target.crystal_model,
                        source.crystal_model,
                    )
                    if zone is not None:
                        predicted_zone[source.key] = zone, error_deg
                        constrain_phase_from_axis(target, source, zone)
                        constrained_key = source.key
            if constrained_key is None:
                warnings.append(
                    "The single direction pair leaves rotation about its common "
                    "axis undetermined, and neither independently fitted phase is "
                    "on the specified axis. Independent image fits were kept."
                )
        else:
            # The paired secondary directions are known non-None here.
            assert orientation.source_secondary is not None
            assert orientation.target_secondary is not None
            rotation = rotation_from_parallel_pairs(
                orientation.source_primary,
                orientation.source_secondary,
                orientation.target_primary,
                orientation.target_secondary,
                source.crystal_structure,
                target.crystal_structure,
                source.crystal_model,
                target.crystal_model,
            )
            if source.key not in selected:
                warnings.append(
                    "Phase 1 did not fit, so its two-pair orientation-relationship "
                    "prediction was not used."
                )
            else:
                source_match = selected[source.key].best
                zone, error_deg = predict_target_zone(
                    rotation,
                    source_match.zone,
                    source.crystal_structure,
                    target.crystal_structure,
                    source.crystal_model,
                    target.crystal_model,
                )
                predicted_zone[target.key] = zone, error_deg
                if error_deg > 2.0:
                    warnings.append(
                        f"The orientation-relationship prediction for {target.label} "
                        f"has a {error_deg:.2f}° residual from the nearest low-index "
                        "direction. Treat it as approximate."
                    )
                constrain_phase_from_axis(source, target, zone)
                fitted_target = selected.get(target.key)
                if (
                    fitted_target is not None
                    and ZAF.canonical_line(fitted_target.best.zone)
                    == ZAF.canonical_line(zone)
                ):
                    projected, predicted_angle = _project_or_fit(
                        source_match,
                        source,
                        fitted_target,
                        target,
                        rotation,
                        rgb.size,
                        alpha_deg,
                        beta_deg,
                        image_to_holder_rotation_deg,
                        holder_order,
                    )
                    angle_error = _line_angle_error_deg(
                        fitted_target.best.angle_deg, predicted_angle
                    )
                    inplane_error[target.key] = angle_error
                    evidence_floor = max(4, int(0.65 * fitted_target.best.matched_count))
                    score_floor = fitted_target.best.score - max(
                        3.0, 0.15 * abs(fitted_target.best.score)
                    )
                    projected_supported = (
                        projected.matched_count >= evidence_floor
                        and projected.score >= score_floor
                    )
                    if angle_error > 5.0 and projected_supported:
                        # The two-pair OR selects the in-plane symmetry variant.
                        ZAF.apply_lattice_calibration(
                            [projected],
                            scale_bar_pixels=scale_bar_pixels,
                            scale_bar_value_inv_nm=scale_bar_value_inv_nm,
                            expected_lattice_parameter_nm=_phase_lattice_nm(target),
                        )
                        selected[target.key] = replace(
                            fitted_target, best=projected, matches=[projected]
                        )
                        inplane_error[target.key] = _line_angle_error_deg(
                            projected.angle_deg, predicted_angle
                        )
                    elif angle_error > 5.0:
                        warnings.append(
                            f"{target.label}: its fitted in-plane spot orientation "
                            f"differs by {angle_error:.1f}° from the specified "
                            "two-pair variant; the image-based fit was kept."
                        )

    phase_results: list[PhaseResult] = []
    for phase in phases:
        indexing = selected.get(phase.key)
        if indexing is None:
            continue
        match = indexing.best
        zone_label = ZAF.format_zone_direction(
            match.zone, match.crystal_structure, match.hcp_four_index
        )
        fitted = ZAF.fitted_diffraction_image(
            rgb,
            match,
            title=f"Fitted {phase.label} {zone_label}",
            highlight_matched_only=True,
        )
        crystal_to_zero_holder, _image_axes, _holder_axes = ZAF.orientation_from_match(
            match,
            alpha_deg,
            beta_deg,
            image_to_holder_rotation_deg,
            holder_order,
        )
        families = ZAF.zone_families_for_structure(
            phase.crystal_structure, phase.crystal_model
        )
        targets = ZAF.target_rows(
            match,
            crystal_to_zero_holder,
            alpha_deg,
            beta_deg,
            families,
            include_opposites=False,
            holder_order=holder_order,
            alpha_limits=alpha_limits,
            beta_limits=beta_limits,
        )
        or_data = predicted_zone.get(phase.key)
        phase_results.append(
            PhaseResult(
                spec=phase,
                indexing=indexing,
                fitted_image=fitted,
                target_rows=targets,
                independent_indexing=independent.get(phase.key),
                or_predicted_zone=or_data[0] if or_data is not None else None,
                or_residual_deg=or_data[1] if or_data is not None else None,
                or_inplane_error_deg=inplane_error.get(phase.key),
            )
        )
    composite = _render_composite(rgb, phase_results)
    return MultiPhaseResult(phase_results, composite, warnings)
