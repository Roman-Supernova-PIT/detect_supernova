import subprocess

from astropy.io import fits
from astropy.table import vstack
from photutils.detection import find_peaks

SOURCE_EXTRACTOR_EXECUTABLE = "source-extractor"
DETECTION_CONFIG = "default.sex"
DETECTION_PARA = "default.param"
DETECTION_FILTER = "default.conv"


def detect(
    difference_path,
    save_path,
    source_extractor_executable=None,
    detection_config=None,
    detection_para=None,
    detection_filter=None,
):
    source_extractor_executable = (
        source_extractor_executable or SOURCE_EXTRACTOR_EXECUTABLE
    )
    detection_config = detection_config or DETECTION_CONFIG
    detection_para = detection_para or DETECTION_PARA
    detection_filter = detection_filter or DETECTION_FILTER
    detection_cmd = [
        source_extractor_executable,
        difference_path,
        "-c",
        detection_config,
        "-PARAMETERS_NAME",
        detection_para,
        "-FILTER_NAME",
        detection_filter,
        "-CATALOG_NAME",
        save_path,
    ]
    result = subprocess.run(detection_cmd, capture_output=True, text=True)

    return result


def score_image_detect(score_image_path, catalog_save_path=None, threshold=10, box_size=11):
    """
    Detect based on the peak pixels in the score image.

    This is equivalent to doing PSF convolution on a direct image and looking for peaks.

    Use astropy.photutils

    Based on
    https://photutils.readthedocs.io/en/stable/user_guide/detection.html
    """
    image = fits.getdata(score_image_path)
    pos_obj = find_peaks(image, threshold=threshold, box_size=box_size)
    neg_obj = find_peaks(-image, threshold=threshold, box_size=box_size)
    neg_obj["peak_value"] = - neg_obj["peak_value"]

    obj = vstack(pos_obj, neg_obj)

    if catalog_save_path is not None:
        obj.write(catalog_save_path)

    return obj
