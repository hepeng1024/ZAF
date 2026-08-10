import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from PIL import Image

import ZAF as finder


PROJECT_DIR = Path(__file__).resolve().parents[1]
KNOWN_110_IMAGE = PROJECT_DIR / "0039_DP_110_zone_840x.bmp"
KNOWN_123_IMAGE = PROJECT_DIR / "BM-OneView_200kV_410X_0071.bmp"
KNOWN_100_IMAGE = PROJECT_DIR / "BM-OneView_200kV_410X_0072.bmp"


class ScaleBarDetectionTests(unittest.TestCase):
    def test_full_width_white_bmp_border_is_not_a_scale_bar(self) -> None:
        pixels = np.zeros((101, 121), dtype=np.uint8)
        pixels[80:84, 10:51] = 255
        pixels[-1, :] = 255

        with TemporaryDirectory() as temporary_directory:
            image_path = Path(temporary_directory) / "bordered.bmp"
            Image.fromarray(pixels).save(image_path)
            detected = finder.detect_scale_bar_pixels(image_path)

        self.assertIsNotNone(detected)
        self.assertEqual(float(detected["length_px"]), 41.0)  # type: ignore[index]


@unittest.skipUnless(
    KNOWN_123_IMAGE.exists() and KNOWN_100_IMAGE.exists(),
    "The supplied 0071/0072 regression images are not available.",
)
class ReportedPatternRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        gray_123, rgb_123 = finder.load_grayscale(KNOWN_123_IMAGE)
        gray_100, rgb_100 = finder.load_grayscale(KNOWN_100_IMAGE)
        common = {
            "current_zone": None,
            "center_xy": None,
            "n_peaks": 120,
            "min_distance_px": None,
            "spot_sigma_px": None,
            "peak_percentile": 99.0,
            "max_index": finder.DEFAULT_REFLECTION_MAX_INDEX,
            "max_g_norm": finder.DEFAULT_REFLECTION_MAX_G_NORM,
            "tolerance_fraction": 0.18,
        }
        cls.indexed_123 = finder.index_diffraction_pattern(
            gray_123,
            rgb_123.size,
            **common,
        )
        cls.indexed_100 = finder.index_diffraction_pattern(
            gray_100,
            rgb_100.size,
            **common,
        )

    def test_0071_uncalibrated_ranking_uses_full_spot_geometry(self) -> None:
        self.assertEqual(self.indexed_123.best.family, "123")
        self.assertAlmostEqual(self.indexed_123.best.scale, 165.917, places=2)
        self.assertLess(self.indexed_123.best.rms_px, 10.0)

    def test_0071_lattice_calibration_keeps_the_primitive_123_fit(self) -> None:
        scale_bar = finder.detect_scale_bar_pixels(KNOWN_123_IMAGE)
        self.assertIsNotNone(scale_bar)

        calibrated = finder.apply_lattice_calibration(
            self.indexed_123.matches,
            scale_bar_pixels=float(scale_bar["length_px"]),  # type: ignore[index]
            scale_bar_value_inv_nm=10.0,
            expected_lattice_parameter_nm=0.415,
        )
        self.assertEqual(calibrated[0].family, "123")
        self.assertAlmostEqual(calibrated[0].scale, 165.917, places=2)
        self.assertAlmostEqual(
            float(calibrated[0].estimated_lattice_nm),
            0.4104,
            places=3,
        )

    def test_0072_uncalibrated_ranking_prefers_100_over_125(self) -> None:
        self.assertEqual(self.indexed_100.best.family, "100")
        families = [match.family for match in self.indexed_100.matches[:2]]
        self.assertEqual(families, ["100", "125"])

    def test_0072_calibration_ignores_export_padding_and_confirms_100(self) -> None:
        scale_bar = finder.detect_scale_bar_pixels(KNOWN_100_IMAGE)
        self.assertIsNotNone(scale_bar)
        self.assertEqual(float(scale_bar["length_px"]), 681.0)  # type: ignore[index]

        calibrated = finder.apply_lattice_calibration(
            self.indexed_100.matches,
            scale_bar_pixels=float(scale_bar["length_px"]),  # type: ignore[index]
            scale_bar_value_inv_nm=10.0,
            expected_lattice_parameter_nm=0.415,
        )
        self.assertEqual(calibrated[0].family, "100")
        self.assertAlmostEqual(
            float(calibrated[0].estimated_lattice_nm),
            0.4115,
            places=3,
        )


@unittest.skipUnless(
    KNOWN_110_IMAGE.exists(),
    "The supplied 110-zone regression image is not available.",
)
class AmbiguityResolutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gray, cls.rgb = finder.load_grayscale(KNOWN_110_IMAGE)

    def test_adaptive_pass_resolves_124_alias_as_110(self) -> None:
        indexing = finder.index_diffraction_pattern(
            self.gray,
            self.rgb.size,
            current_zone=None,
            center_xy=None,
            n_peaks=120,
            min_distance_px=None,
            spot_sigma_px=None,
            peak_percentile=99.0,
            max_index=finder.DEFAULT_REFLECTION_MAX_INDEX,
            max_g_norm=finder.DEFAULT_REFLECTION_MAX_G_NORM,
            tolerance_fraction=0.18,
        )

        self.assertEqual(indexing.diagnostics.initial_best_family, "124")
        self.assertTrue(indexing.diagnostics.adaptive_used)
        self.assertEqual(indexing.diagnostics.initial_peak_count, 24)
        self.assertGreater(indexing.diagnostics.final_peak_count, 24)
        self.assertEqual(indexing.best.family, "110")
        self.assertFalse(indexing.diagnostics.ambiguous)

    def test_scale_and_known_lattice_resolve_sparse_match(self) -> None:
        peaks = finder.detect_spots(
            self.gray,
            n_peaks=120,
            peak_percentile=99.0,
        )
        best, matches = finder.choose_best_match(
            peaks,
            self.rgb.size,
            current_zone=None,
            center_xy=None,
            max_index=finder.DEFAULT_REFLECTION_MAX_INDEX,
            max_g_norm=finder.DEFAULT_REFLECTION_MAX_G_NORM,
            tolerance_fraction=0.18,
        )
        ambiguous, _gap, _required_gap = finder.match_ambiguity(matches)
        self.assertEqual(best.family, "124")
        self.assertTrue(ambiguous)

        scale_bar = finder.detect_scale_bar_pixels(KNOWN_110_IMAGE)
        self.assertIsNotNone(scale_bar)
        calibrated = finder.apply_lattice_calibration(
            matches,
            scale_bar_pixels=float(scale_bar["length_px"]),  # type: ignore[index]
            scale_bar_value_inv_nm=5.0,
            expected_lattice_parameter_nm=0.3524,
        )

        self.assertEqual(calibrated[0].family, "110")
        self.assertAlmostEqual(
            float(calibrated[0].estimated_lattice_nm),
            0.3487,
            places=3,
        )


if __name__ == "__main__":
    unittest.main()
