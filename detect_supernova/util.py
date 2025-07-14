from dataclasses import dataclass
import os

import pathlib


INPUT_IMAGE_PATTERN = (
    "RomanTDS/images/simple_model/{band}/{pointing}/Roman_TDS_simple_model_{band}_{pointing}_{sca}.fits.gz"
)
INPUT_TRUTH_PATTERN = (
    "RomanTDS/truth/{band}/{pointing}/Roman_TDS_index_{band}_{pointing}_{sca}.txt"
)
SIMS_DIR = pathlib.Path(os.getenv("SIMS_DIR", None))


IMAGE_WIDTH = 4088
IMAGE_HEIGHT = 4088


@dataclass
class ImageInfo:
    data_id: dict
    temp_dir: pathlib.Path

    def __post_init__(self):
        self.image_path = SIMS_DIR / pathlib.Path(
            INPUT_IMAGE_PATTERN.format(**self.data_id)
        )
        self.cx = IMAGE_WIDTH // 2
        self.cy = IMAGE_HEIGHT // 2

        self.image_name = self.image_path.name
        self.skysub_path = self.temp_dir / f"skysub_{self.image_name}"
        self.detmask_path = self.temp_dir / f"detmask_{self.image_name}"
        self.psf_path = self.temp_dir / f"psf_{self.image_name}"
