import unittest

import numpy as np

import ZAF as finder
import ZAF_gui as gui


class CrystalStructureTests(unittest.TestCase):
    def test_structure_specific_reflection_conditions(self) -> None:
        self.assertTrue(finder.reflection_allowed("FCC", 1, 1, 1))
        self.assertFalse(finder.reflection_allowed("FCC", 1, 0, 0))

        self.assertTrue(finder.reflection_allowed("BCC", 1, 1, 0))
        self.assertFalse(finder.reflection_allowed("BCC", 1, 0, 0))

        self.assertTrue(finder.reflection_allowed("HCP", 0, 0, 2))
        self.assertFalse(finder.reflection_allowed("HCP", 0, 0, 1))
        self.assertTrue(finder.reflection_allowed("HCP", 1, 0, 0))

    def test_hcp_direct_and_reciprocal_bases_are_dual(self) -> None:
        direction = (1, 2, 3)
        reflection = (-2, 4, 1)
        direct = finder.direction_cartesian(direction, "HCP", 1.58)
        reciprocal = finder.reciprocal_cartesian(reflection, "HCP", 1.58)
        self.assertAlmostEqual(
            float(np.dot(direct, reciprocal)),
            float(np.dot(direction, reflection)),
            places=12,
        )

    def test_hcp_three_and_four_index_direction_round_trip(self) -> None:
        for direction in ((1, 0, 0), (1, 2, 3), (-2, 1, 4), (0, 0, 1)):
            with self.subTest(direction=direction):
                four_index = finder.three_to_four_direction(direction)
                self.assertEqual(sum(four_index[:3]), 0)
                self.assertEqual(
                    finder.four_to_three_direction(four_index),
                    finder.reduce_miller(direction),
                )

        self.assertEqual(
            finder.parse_zone_direction("[2 -1 -1 0]", "HCP", True),
            (1, 0, 0),
        )
        self.assertEqual(
            finder.parse_zone_direction("2-1-10", "HCP", True),
            (1, 0, 0),
        )
        self.assertEqual(
            finder.format_zone_direction((1, 0, 0), "HCP", True),
            "[2-1-10]",
        )

    def test_hcp_uses_the_dedicated_nonoverlapping_family_catalog(self) -> None:
        expected = (
            "0001",
            "2-1-10",
            "10-10",
            "10-11",
            "10-12",
            "11-23",
            "21-30",
            "40-43",
        )
        self.assertEqual(finder.HCP_ZONE_FAMILIES, expected)
        self.assertEqual(finder.zone_families_for_structure("HCP"), expected)
        self.assertTrue(set(expected).issubset(finder.ZONE_FAMILY_COLORS))

        seen: set[tuple[int, int, int]] = set()
        for family in expected:
            directions = set(
                finder.family_directions(family, crystal_structure="HCP")
            )
            self.assertTrue(directions)
            self.assertFalse(seen.intersection(directions))
            self.assertTrue(
                all(
                    finder.family_name(direction, "HCP") == family
                    for direction in directions
                )
            )
            seen.update(directions)

    def test_tilt_simulator_hcp_cell_matches_crysdis_basis(self) -> None:
        atoms, edge_points = gui.SampleTiltSimulator.create_crysdis_hcp_lattice(
            size=5.6,
            hcp_c_over_a=finder.DEFAULT_HCP_C_OVER_A,
        )
        self.assertEqual(atoms.shape, (3, 9))
        self.assertEqual(edge_points.shape, (3, 24))
        self.assertEqual(
            len({tuple(np.round(point, 8)) for point in atoms.T}),
            9,
        )

    def test_generated_reflections_obey_the_zone_law(self) -> None:
        zone = (1, 2, 3)
        for structure in finder.CRYSTAL_STRUCTURES:
            with self.subTest(structure=structure):
                reflections, _basis_x, _basis_y = finder.make_reflections(
                    zone,
                    structure,
                    max_index=6,
                    max_g_norm=8.0,
                )
                self.assertGreater(len(reflections), 1)
                for reflection in reflections:
                    self.assertEqual(
                        reflection.h * zone[0]
                        + reflection.k * zone[1]
                        + reflection.l * zone[2],
                        0,
                    )
                    self.assertTrue(
                        reflection.h == reflection.k == reflection.l == 0
                        or finder.reflection_allowed(
                            structure,
                            reflection.h,
                            reflection.k,
                            reflection.l,
                        )
                    )

    def test_ideal_bcc_and_hcp_patterns_auto_detect(self) -> None:
        image_size = (768, 768)
        center = (384.0, 384.0)
        cases = (("BCC", "123"), ("HCP", "10-11"))
        for structure, family in cases:
            with self.subTest(structure=structure, family=family):
                reflections, _basis_x, _basis_y = finder.make_reflections(
                    finder.parse_family_direction(family, structure),
                    structure,
                    max_index=8,
                    max_g_norm=9.0,
                )
                screen_points = (
                    np.asarray([reflection.xy for reflection in reflections]) * 38.0
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
                    max_index=8,
                    max_g_norm=9.0,
                    tolerance_fraction=0.18,
                    crystal_structure=structure,
                )
                self.assertEqual(best.crystal_structure, structure)
                self.assertEqual(best.family, family)

    def test_ideal_hcp_patterns_auto_detect_each_dedicated_family(self) -> None:
        image_size = (1024, 1024)
        center = (512.0, 512.0)
        for family in finder.HCP_ZONE_FAMILIES:
            with self.subTest(family=family):
                zone = finder.parse_family_direction(family, "HCP")
                reflections, _basis_x, _basis_y = finder.make_reflections(
                    zone,
                    "HCP",
                    max_index=10,
                    max_g_norm=11.0,
                )
                screen_points = (
                    np.asarray([reflection.xy for reflection in reflections]) * 36.0
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
                    max_index=10,
                    max_g_norm=11.0,
                    tolerance_fraction=0.18,
                    crystal_structure="HCP",
                    hcp_four_index=True,
                )
                self.assertEqual(best.family, family)


if __name__ == "__main__":
    unittest.main()
