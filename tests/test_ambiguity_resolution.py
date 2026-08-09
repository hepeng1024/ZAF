import unittest
from pathlib import Path

import ZAF as finder


PROJECT_DIR = Path(__file__).resolve().parents[1]
KNOWN_110_IMAGE = PROJECT_DIR / "0039_DP_110_zone_840x.bmp"


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
