import subprocess

SOURCE_EXTRACTOR_EXECUTABLE = "sex"
DETECTION_CONFIG = "default.sex"
DETECTION_PARA = "default.param"
DETECTION_FILTER = "default.conv"

def detect(difference_path, save_path, source_extractor_executable=None, detection_config=None, detection_para=None, detection_filter=None):
    source_extractor_executable = source_extractor_executable or SOURCE_EXTRACTOR_EXECUTABLE
    detection_config = detection_config or DETECTION_CONFIG
    detection_para = detection_para or DETECTION_PARA
    detection_filter = detection_filter or DETECTION_FILTER
    detection_cmd = [source_extractor_executable, difference_path,
                     "-c", detection_config, "-PARAMETERS_NAME", detection_para,
                     "-FILTER_NAME", detection_filter, "-CATALOG_NAME", save_path]
    result = subprocess.run(detection_cmd, capture_output=True, text=True)