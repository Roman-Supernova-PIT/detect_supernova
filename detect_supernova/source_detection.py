import subprocess

from astropy.io import fits
from astropy.table import vstack
from photutils.centroids import centroid_com
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


def score_image_detect(
    score_image_path, catalog_save_path=None, threshold=10, box_size=11, negative=True, overwrite=True,
):
    """Detect based on the peak pixels in the score image.

    Parameters
    ----------
    score_image_path : str
        Path to score image
    threshold : float
        Signal-to-noise ratio threshold.
    box_size : int
        Size of box in which to look for unique peaks.  Passed to photutils.find_peaks.
    negative : bool
        Search for negative sources as well as positive sources.
    overwrite : bool
        Overwrite existing catalog_save_path

    Returns
    -------
    AstroPy Table of results

    Notes
    -----
    The score image for a subtraction is the equivalent of cross-correlating a
    direct image with its PSF and dividing by the variance.  We then look for
    for significant peaks to identiy sources.

    The searches for positive and negative sources run separately,
    and thus a positive source and a negative source can be found within the
    same box_size region.

    Uses astropy.photutils

    Based on
    https://photutils.readthedocs.io/en/stable/user_guide/detection.html
    """
    image = fits.getdata(score_image_path)
    find_peaks_kwargs = {"threshold": threshold, "box_size": box_size, "centroid_func": centroid_com}
    pos_obj = find_peaks(image, **find_peaks_kwargs)
    neg_obj = find_peaks(-image, **find_peaks_kwargs)
    neg_obj["peak_value"] = -neg_obj["peak_value"]

    obj = vstack([pos_obj, neg_obj])

    if catalog_save_path is not None:
        obj.write(catalog_save_path, overwrite=overwrite)

    return obj
