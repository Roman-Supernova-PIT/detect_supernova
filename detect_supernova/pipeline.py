import argparse
import atexit
import os
import pathlib
import tempfile

import pandas as pd

from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_skycoord

import data_loader
import subtraction
import source_detection
import truth_matching
import truth_retrieval


class Detection:

    INPUT_COLUMNS = [
        "science_band",
        "science_pointing",
        "science_sca",
        "template_band",
        "template_pointing",
        "template_sca",
    ]

    """
    INPUT_IMAGE_PATTERN = ("/global/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data"
                                "/RomanTDS/images/simple_model/{band}/{pointing}/Roman_TDS_simple_model_{band}_{pointing}_{sca}.fits.gz")
    INPUT_TRUTH_PATTERN = ("/global/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data"
                                 "/RomanTDS/truth/{band}/{pointing}/Roman_TDS_index_{band}_{pointing}_{sca}.txt")
    """
    SIMS_DIR = os.getenv("SIMS_DIR", None)

    INPUT_IMAGE_PATTERN = (
        SIMS_DIR
        + "/RomanTDS/images/simple_model/{band}/{pointing}/Roman_TDS_simple_model_{band}_{pointing}_{sca}.fits.gz"
    )
    INPUT_TRUTH_PATTERN = (
        SIMS_DIR
        + "/RomanTDS/truth/{band}/{pointing}/Roman_TDS_index_{band}_{pointing}_{sca}.txt"
    )

    DIFF_PATTERN = "{science_band}_{science_pointing}_{science_sca}_-_{template_band}_{template_pointing}_{template_sca}"

    # Source detection config.
    SOURCE_EXTRACTOR_EXECUTABLE = "source-extractor"
    DETECTION_CONFIG = os.path.join(
        os.path.dirname(__file__), "..", "configs", "default.sex"
    )
    DETECTION_PARA = os.path.join(
        os.path.dirname(__file__), "..", "configs", "default.param"
    )
    DETECTION_FILTER = os.path.join(
        os.path.dirname(__file__), "..", "configs", "default.conv"
    )

    # Source Matching
    MATCH_RADIUS = 0.4  # in arcsec unit

    # file prefix
    DIFF_IMAGE_PREFIX = "decorr_diff_"
    DIFF_SCORE_PREFIX = "score_"
    DIFF_DETECTION_PREFIX = "detection_"
    DIFF_SCORE_DETECTION_PREFIX = "score_detection_"
    DIFF_TRUTH_PREFIX = "truth_"
    TRANSIENTS_TO_DETECTION_PREFIX = "transients_to_detection_"
    DETECTION_TO_TRANSIENTS_PREFIX = "detection_to_transients_"
    TRANSIENTS_TO_SCORE_DETECTION_PREFIX = "transients_to_score_detection_"
    SCORE_DETECTION_TO_TRANSIENTS_PREFIX = "score_detection_to_transients_"

    def __init__(self, data_records_path, temp_dir=None, output_dir="./output"):
        self.data_records_path = data_records_path
        self.data_records = pd.read_csv(
            self.data_records_path, usecols=self.INPUT_COLUMNS
        )
        self.temp_dir = temp_dir
        self.output_dir = output_dir

    @staticmethod
    def retrieve_truth(
        science_image_path,
        template_image_path,
        science_truth_path,
        template_truth_path,
        difference_truth_path,
    ):
        science_wcs = data_loader.load_wcs(science_image_path, hdu_id=1)
        template_wcs = data_loader.load_wcs(template_image_path, hdu_id=1)

        science_truth = data_loader.load_table(science_truth_path)
        template_truth = data_loader.load_table(template_truth_path)

        truth = truth_retrieval.merge_science_and_template_truth(
            science_truth, template_truth, science_wcs, template_wcs, offset=50
        )
        truth.to_csv(difference_truth_path, index=False)
        return truth

    @staticmethod
    def match_transients(
        truth,
        difference_image_path,
        difference_detection_path,
        match_radius,
        transients_to_detection_path,
        detection_to_transients_path,
    ):
        difference_wcs = data_loader.load_wcs(difference_image_path, hdu_id=0)

        detection = data_loader.load_table(difference_detection_path)
        transients = truth[truth.obj_type == "transient"].copy().reset_index(drop=True)
        transients_skycoord = SkyCoord(
            transients.ra, transients.dec, frame="icrs", unit="deg"
        )
        detection_skycoord = pixel_to_skycoord(
            detection.X_IMAGE, detection.Y_IMAGE, difference_wcs
        )
        transients_to_detection = truth_matching.skymatch_and_join(
            transients, detection, transients_skycoord, detection_skycoord, match_radius
        )
        detection_to_transients = truth_matching.skymatch_and_join(
            detection, transients, detection_skycoord, transients_skycoord, match_radius
        )

        transients_to_detection.to_csv(transients_to_detection_path, index=False)
        detection_to_transients.to_csv(detection_to_transients_path, index=False)
        return transients_to_detection, detection_to_transients

    @staticmethod
    def get_difference_id(science_id, template_id):
        _prefixed_science = {f"science_{k}": v for k, v in science_id.items()}
        _prefixed_template = {f"template_{k}": v for k, v in template_id.items()}
        difference_id = {**_prefixed_science, **_prefixed_template}
        return difference_id

    def path_helper(self, science_id, template_id):
        file_path = {}

        difference_id = self.__class__.get_difference_id(science_id, template_id)
        diff_pattern = self.DIFF_PATTERN.format(**difference_id)
        file_path["full_output_dir"] = os.path.join(self.output_dir, diff_pattern)
        os.makedirs(file_path["full_output_dir"], exist_ok=True)

        # subtraction
        file_path["science_image_path"] = self.INPUT_IMAGE_PATTERN.format(**science_id)
        file_path["template_image_path"] = self.INPUT_IMAGE_PATTERN.format(
            **template_id
        )
        file_path["difference_image_path"] = os.path.join(
            file_path["full_output_dir"],
            self.DIFF_IMAGE_PREFIX + diff_pattern + ".fits",
        )
        # detection
        file_path["difference_detection_path"] = os.path.join(
            file_path["full_output_dir"],
            self.DIFF_DETECTION_PREFIX + diff_pattern + ".cat",
        )
        file_path["score_image_path"] = os.path.join(
            file_path["full_output_dir"],
            self.DIFF_SCORE_PREFIX + diff_pattern + ".fits",
        )
        # diff score image detection
        file_path["score_image_detection_path"] = os.path.join(
            file_path["full_output_dir"],
            self.SCORE_DETECTION_PREFIX + diff_pattern + ".ecsv",
        )
        # decorr diff image detection
        file_path["difference_detection_path"] = os.path.join(
            file_path["full_output_dir"],
            self.DIFF_DETECTION_PREFIX + diff_pattern + ".cat",
        )
        # truth retrieval
        file_path["science_truth_path"] = self.INPUT_TRUTH_PATTERN.format(**science_id)
        file_path["template_truth_path"] = self.INPUT_TRUTH_PATTERN.format(
            **template_id
        )
        file_path["difference_truth_path"] = os.path.join(
            file_path["full_output_dir"],
            self.DIFF_TRUTH_PREFIX + diff_pattern + ".fits",
        )
        # truth matching
        file_path["transients_to_detection_path"] = os.path.join(
            file_path["full_output_dir"],
            self.TRANSIENTS_TO_DETECTION_PREFIX + diff_pattern + ".csv",
        )
        file_path["detection_to_transients_path"] = os.path.join(
            file_path["full_output_dir"],
            self.DETECTION_TO_TRANSIENTS_PREFIX + diff_pattern + ".csv",
        )
        file_path["transients_to_score_detection_path"] = os.path.join(
            file_path["full_output_dir"],
            self.TRANSIENTS_TO_SCORE_DETECTION_PREFIX + diff_pattern + ".csv",
        )
        file_path["score_detection_to_transients_path"] = os.path.join(
            file_path["full_output_dir"],
            self.SCORE_DETECTION_TO_TRANSIENTS_PREFIX + diff_pattern + ".csv",
        )
        return file_path

    def run_one_subtraction(
        self,
        science_band,
        science_pointing,
        science_sca,
        template_band,
        template_pointing,
        template_sca,
        temp_dir,
    ):
        science_id = {
            "band": science_band,
            "pointing": science_pointing,
            "sca": science_sca,
        }
        template_id = {
            "band": template_band,
            "pointing": template_pointing,
            "sca": template_sca,
        }
        file_path = self.path_helper(science_id, template_id)

        print(
            "[INFO] Processing started for data records "
            f"| Science ID {science_id} "
            f"| Template ID {template_id} "
        )

        print("[INFO] Processing subtraction")
        subtract = subtraction.Pipeline(
            science_band=science_band,
            science_pointing=science_pointing,
            science_sca=science_sca,
            template_band=template_band,
            template_pointing=template_pointing,
            template_sca=template_sca,
            temp_dir=temp_dir,
            out_dir=file_path["full_output_dir"],
        )
        subtract.run()

        print("[INFO] Processing detection")
        source_detection.detect(
            file_path["difference_image_path"],
            file_path["difference_detection_path"],
            source_extractor_executable=self.SOURCE_EXTRACTOR_EXECUTABLE,
            detection_config=self.DETECTION_CONFIG,
            detection_para=self.DETECTION_PARA,
            detection_filter=self.DETECTION_FILTER,
        )

        print("[INFO] Processing score image detection")
        source_detection.score_image_detect(
            file_path["score_image_path"],
            file_path["score_image_detection_path"],
        )

        print("[INFO] Processing truth retrieval")
        truth = self.__class__.retrieve_truth(
            file_path["science_image_path"],
            file_path["template_image_path"],
            file_path["science_truth_path"],
            file_path["template_truth_path"],
            file_path["difference_truth_path"],
        )

        print("[INFO] Processing diffim detection truth matching")
        _, _ = self.__class__.match_transients(
            truth,
            file_path["difference_image_path"],
            file_path["difference_detection_path"],
            self.MATCH_RADIUS,
            file_path["transients_to_detection_path"],
            file_path["detection_to_transients_path"],
        )

        print("[INFO] Processing score image detection truth matching")
        _, _ = self.__class__.match_transients(
            truth,
            file_path["difference_image_path"],
            file_path["score_image_detection_path"],
            self.MATCH_RADIUS,
            file_path["transients_to_score_image_detection_path"],
            file_path["score_image_detection_to_transients_path"],
        )

        print("[INFO] Processing finished.")

    def run(self):
        os.makedirs(self.output_dir, exist_ok=True)

        # create temporary directory
        if self.temp_dir is None:
            temp_dir_obj = tempfile.TemporaryDirectory()
            temp_dir = pathlib.Path(temp_dir_obj.name)
            atexit.register(temp_dir_obj.cleanup)
        else:
            temp_dir = pathlib.Path(self.temp_dir)
            os.makedirs(temp_dir, exist_ok=True)

        for i, row in self.data_records.iterrows():
            self.run_one_subtraction(
                row["science_band"],
                row["science_pointing"],
                row["science_sca"],
                row["template_band"],
                row["template_pointing"],
                row["template_sca"],
                temp_dir=temp_dir,
            )


def main():
    parser = argparse.ArgumentParser("detection pipeline")
    parser.add_argument(
        "-d",
        "--data-records",
        type=str,
        required=True,
        help="Input file with data records.",
    )
    parser.add_argument("-t", "--temp-dir", default=None, help="Temporary directory.")
    parser.add_argument(
        "-o", "--output-dir", type=str, default="./output", help="Output directory."
    )
    args = parser.parse_args()

    detection = Detection(args.data_records, args.temp_dir, args.output_dir)
    detection.run()


if __name__ == "__main__":
    main()
