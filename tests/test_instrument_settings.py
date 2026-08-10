from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import ZAF_gui


VALID_SETTINGS = """
# Institution-specific TEM values
alpha_min = -40
alpha_max = 38.5
beta_min = -25
beta_max = 22
image_to_holder_rotation_deg = 87.25  # camera calibration
"""


class InstrumentSettingsTests(unittest.TestCase):
    def test_repository_template_is_valid(self) -> None:
        parsed = ZAF_gui.parse_instrument_settings(
            ZAF_gui.BUNDLED_INSTRUMENT_SETTINGS_PATH.read_text(encoding="utf-8")
        )
        self.assertEqual(parsed, ZAF_gui.DEFAULT_INSTRUMENT_SETTINGS)

    def test_parser_accepts_comments_whitespace_and_decimal_values(self) -> None:
        parsed = ZAF_gui.parse_instrument_settings(VALID_SETTINGS)
        self.assertEqual(parsed["alpha_min"], -40.0)
        self.assertEqual(parsed["alpha_max"], 38.5)
        self.assertEqual(parsed["beta_min"], -25.0)
        self.assertEqual(parsed["beta_max"], 22.0)
        self.assertEqual(parsed["image_to_holder_rotation_deg"], 87.25)

    def test_parser_rejects_invalid_files(self) -> None:
        invalid_cases = {
            "missing key": VALID_SETTINGS.replace("beta_max = 22\n", ""),
            "unknown key": VALID_SETTINGS + "holder_order = 1\n",
            "duplicate key": VALID_SETTINGS + "alpha_min = -30\n",
            "nonnumeric": VALID_SETTINGS.replace("alpha_min = -40", "alpha_min = left"),
            "nonfinite": VALID_SETTINGS.replace("alpha_min = -40", "alpha_min = nan"),
            "reversed alpha": VALID_SETTINGS.replace("alpha_min = -40", "alpha_min = 50"),
            "reversed beta": VALID_SETTINGS.replace("beta_min = -25", "beta_min = 30"),
            "missing equals": VALID_SETTINGS + "not an assignment\n",
        }
        for label, text in invalid_cases.items():
            with self.subTest(label=label):
                with self.assertRaises(ValueError):
                    ZAF_gui.parse_instrument_settings(text)

    def test_loader_uses_built_in_defaults_and_warning_on_error(self) -> None:
        with tempfile.TemporaryDirectory() as temp_name:
            temp_dir = Path(temp_name)
            missing_values, missing_warning = ZAF_gui.load_instrument_settings(
                temp_dir / "missing.txt"
            )
            self.assertEqual(
                missing_values, ZAF_gui.DEFAULT_INSTRUMENT_SETTINGS
            )
            self.assertIn("missing", missing_warning or "")

            invalid_path = temp_dir / "invalid.txt"
            invalid_path.write_text("alpha_min = broken\n", encoding="utf-8")
            invalid_values, invalid_warning = ZAF_gui.load_instrument_settings(
                invalid_path
            )
            self.assertEqual(
                invalid_values, ZAF_gui.DEFAULT_INSTRUMENT_SETTINGS
            )
            self.assertIn("invalid", invalid_warning or "")

    def test_loader_returns_valid_custom_values_without_warning(self) -> None:
        with tempfile.TemporaryDirectory() as temp_name:
            path = Path(temp_name) / "settings.txt"
            path.write_text(VALID_SETTINGS, encoding="utf-8")
            values, warning = ZAF_gui.load_instrument_settings(path)
        self.assertEqual(values["image_to_holder_rotation_deg"], 87.25)
        self.assertIsNone(warning)

    def test_frozen_linux_settings_file_is_beside_executable(self) -> None:
        executable = "/opt/ZAF release/ZAF"
        with (
            patch.object(ZAF_gui.sys, "frozen", True, create=True),
            patch.object(ZAF_gui.sys, "platform", "linux"),
            patch.object(ZAF_gui.sys, "executable", executable),
        ):
            path = ZAF_gui.instrument_settings_path()
        self.assertEqual(
            path, Path(executable).parent / ZAF_gui.INSTRUMENT_SETTINGS_FILENAME
        )

    def test_frozen_macos_settings_file_uses_application_support(self) -> None:
        executable = (
            "/private/var/folders/example/AppTranslocation/ABC/d/"
            "ZAF.app/Contents/MacOS/ZAF"
        )
        with (
            patch.object(ZAF_gui.sys, "frozen", True, create=True),
            patch.object(ZAF_gui.sys, "platform", "darwin"),
            patch.object(ZAF_gui.sys, "executable", executable),
        ):
            path = ZAF_gui.instrument_settings_path(Path("/Users/nina"))
        self.assertEqual(
            path,
            Path("/Users/nina/Library/Application Support/ZAF")
            / ZAF_gui.INSTRUMENT_SETTINGS_FILENAME,
        )

    def test_macos_initialization_copies_bundled_template_once(self) -> None:
        with tempfile.TemporaryDirectory() as temp_name:
            temp_dir = Path(temp_name)
            template = temp_dir / "bundle" / ZAF_gui.INSTRUMENT_SETTINGS_FILENAME
            settings = (
                temp_dir
                / "Library"
                / "Application Support"
                / "ZAF"
                / ZAF_gui.INSTRUMENT_SETTINGS_FILENAME
            )
            template.parent.mkdir()
            template.write_text(VALID_SETTINGS, encoding="utf-8")

            warning = ZAF_gui.initialize_instrument_settings(settings, template)
            self.assertIsNone(warning)
            self.assertEqual(settings.read_text(encoding="utf-8"), VALID_SETTINGS)

            custom_settings = VALID_SETTINGS.replace(
                "alpha_min = -40", "alpha_min = -45"
            )
            settings.write_text(custom_settings, encoding="utf-8")
            warning = ZAF_gui.initialize_instrument_settings(settings, template)
            self.assertIsNone(warning)
            self.assertEqual(settings.read_text(encoding="utf-8"), custom_settings)


if __name__ == "__main__":
    unittest.main()
