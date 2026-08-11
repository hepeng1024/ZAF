import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib
from PIL import ImageFont

import ZAF


class RenderFontTests(unittest.TestCase):
    def test_matplotlib_bold_font_is_the_cross_platform_first_choice(self) -> None:
        expected = (
            Path(matplotlib.get_data_path())
            / "fonts"
            / "ttf"
            / "DejaVuSans-Bold.ttf"
        )

        self.assertTrue(expected.is_file())
        self.assertEqual(ZAF._label_font_candidates()[0], expected)

    def test_label_font_preserves_requested_pixel_size(self) -> None:
        font = ZAF.load_label_font(72)
        bbox = font.getbbox("123")

        self.assertIsInstance(font, ImageFont.FreeTypeFont)
        self.assertGreater(bbox[3] - bbox[1], 35)

    def test_embedded_fallback_is_scaled_when_no_font_file_exists(self) -> None:
        with patch.object(ZAF, "_label_font_candidates", return_value=()):
            font = ZAF.load_label_font(72)
        bbox = font.getbbox("123")

        self.assertGreater(bbox[3] - bbox[1], 35)


if __name__ == "__main__":
    unittest.main()
