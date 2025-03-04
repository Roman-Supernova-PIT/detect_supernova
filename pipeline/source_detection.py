import subprocess

SOURCE_EXTRACTOR_EXECUTABLE = "sex"
DETECTION_CONFIG = "default.sex"

def detect(difference_path, save_path, source_extractor_executable=None, detection_config=None):
    source_extractor_executable = source_extractor_executable or SOURCE_EXTRACTOR_EXECUTABLE
    detection_config = detection_config or DETECTION_CONFIG
    detection_cmd = [source_extractor_executable, difference_path, "-c", "default.sex", "-CATALOG_NAME", save_path]
    result = subprocess.run(detection_cmd, capture_output=True, text=True)