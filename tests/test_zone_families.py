import itertools
import math
import unittest

import numpy as np

import ZAF as finder
import ZAF_gui as gui


ADDED_ZONE_FAMILIES = (
    "105",
    "106",
    "107",
    "115",
    "116",
    "122",
    "123",
    "124",
    "125",
    "133",
    "134",
    "203",
    "205",
    "223",
    "233",
    "304",
    "305",
)


class ZoneFamilyTests(unittest.TestCase):
    def test_supported_families_exhaust_index_sums_through_eight(self) -> None:
        expected_bases = {
            indices
            for indices in itertools.combinations_with_replacement(range(9), 3)
            if 0 < sum(indices) <= 8 and math.gcd(*indices) == 1
        }
        actual_bases = {
            tuple(sorted(finder.parse_family_base(family)))
            for family in finder.SUPPORTED_ZONE_FAMILIES
        }
        self.assertEqual(actual_bases, expected_bases)
        self.assertEqual(len(finder.SUPPORTED_ZONE_FAMILIES), len(expected_bases))
        self.assertTrue(
            set(finder.SUPPORTED_ZONE_FAMILIES).issubset(
                finder.ZONE_FAMILY_COLORS
            )
        )

    def test_added_families_are_registered_with_colors_and_unchecked(self) -> None:
        for family in ADDED_ZONE_FAMILIES:
            with self.subTest(family=family):
                self.assertIn(family, finder.SUPPORTED_ZONE_FAMILIES)
                self.assertIn(family, finder.ZONE_FAMILY_COLORS)
                self.assertNotIn(family, gui.DEFAULT_SELECTED_TARGET_FAMILIES)

    def test_added_series_generate_all_cubic_symmetry_equivalents(self) -> None:
        expected_line_counts = {
            "105": 12,
            "106": 12,
            "107": 12,
            "115": 12,
            "116": 12,
            "122": 12,
            "123": 24,
            "124": 24,
            "125": 24,
            "133": 12,
            "134": 24,
            "203": 12,
            "205": 12,
            "223": 12,
            "233": 12,
            "304": 12,
            "305": 12,
        }
        for family, expected_count in expected_line_counts.items():
            with self.subTest(family=family):
                directions = finder.family_directions(family)
                self.assertEqual(len(directions), expected_count)
                self.assertTrue(all(finder.family_name(direction) == family for direction in directions))
                self.assertEqual(
                    len(finder.family_directions(family, include_opposites=True)),
                    2 * expected_count,
                )

    def test_ideal_patterns_auto_detect_each_added_family(self) -> None:
        image_size = (1024, 1024)
        center = (512.0, 512.0)
        scale = 30.0

        for family in ADDED_ZONE_FAMILIES:
            with self.subTest(family=family):
                reflections, _basis_x, _basis_y = finder.make_fcc_reflections(
                    finder.parse_family_base(family),
                    max_index=finder.DEFAULT_REFLECTION_MAX_INDEX,
                    max_g_norm=finder.DEFAULT_REFLECTION_MAX_G_NORM,
                )
                screen_points = (
                    np.asarray([reflection.xy for reflection in reflections]) * scale
                    + np.asarray([center[0], -center[1]])
                )
                visible = finder.visible_mask(screen_points, *image_size)
                peaks = [
                    finder.Peak(float(point[0]), float(-point[1]), 1.0)
                    for point in screen_points[visible]
                ]

                best, _matches = finder.choose_best_match(
                    peaks,
                    image_size,
                    current_zone=None,
                    center_xy=center,
                    max_index=finder.DEFAULT_REFLECTION_MAX_INDEX,
                    max_g_norm=finder.DEFAULT_REFLECTION_MAX_G_NORM,
                    tolerance_fraction=0.18,
                )

                self.assertEqual(best.family, family)


if __name__ == "__main__":
    unittest.main()
