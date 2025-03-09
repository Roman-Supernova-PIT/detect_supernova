import argparse
import os
import pathlib
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_skycoord
import truth_retrieval
import data_loader
import coord_projection
import subtraction
import source_detection
import truth_matching

class Detection:

    INPUT_COLUMNS = ['science_band', 'science_pointing', 'science_sca', 'template_band', 'template_pointing', 'template_sca']
    
    """
    INPUT_IMAGE_PATTERN = ("/global/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data"
                                "/RomanTDS/images/simple_model/{band}/{pointing}/Roman_TDS_simple_model_{band}_{pointing}_{sca}.fits.gz")
    INPUT_TRUTH_PATTERN = ("/global/cfs/cdirs/lsst/shared/external/roman-desc-sims/Roman_data"
                                 "/RomanTDS/truth/{band}/{pointing}/Roman_TDS_index_{band}_{pointing}_{sca}.txt")
    """
    SIMS_DIR = os.getenv( 'SIMS_DIR', None )
    
    INPUT_IMAGE_PATTERN = SIMS_DIR + "/RomanTDS/images/simple_model/{band}/{pointing}/Roman_TDS_simple_model_{band}_{pointing}_{sca}.fits.gz"
    INPUT_TRUTH_PATTERN = SIMS_DIR + "/RomanTDS/truth/{band}/{pointing}/Roman_TDS_index_{band}_{pointing}_{sca}.txt"

    SUB_OUTPUT_DIR_PATTERN = 'diff_{science_band}_{science_pointing}_{science_sca}_-_{template_band}_{template_pointing}_{template_sca}'

    # Source detection config.
    SOURCE_EXTRACTOR_EXECUTABLE = "sex"
    DETECTION_CONFIG = "default.sex"

    MATCH_RADIUS = 0.4 # in arcsec unit
    
    def __init__(self, data_records_path, base_output_dir='./output'):
        self.data_records_path = data_records_path
        self.data_records = pd.read_csv(self.data_records_path, usecols=self.INPUT_COLUMNS)
        self.base_output_dir = base_output_dir

    @staticmethod
    def get_difference_id(science_id, template_id):
        _prefixed_science = {f"science_{k}": v for k, v in science_id.items()}
        _prefixed_template = {f"template_{k}": v for k, v in template_id.items()}
        difference_id = {**_prefixed_science, **_prefixed_template}
        return difference_id
        
    @staticmethod
    def retrieve_truth(science_image_path, template_image_path,
                       science_truth_path, template_truth_path,
                       difference_truth_path):
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
    def match_transients(truth, difference_image_path, difference_detection_path, match_radius,
                      transients_to_detection_path, detection_to_transients_path):
        difference_wcs = data_loader.load_wcs(difference_image_path, hdu_id=0)
         
        detection = data_loader.load_table(difference_detection_path)
        transients = truth[truth.obj_type=='transient'].copy().reset_index(drop=True)
        transients_skycoord =  SkyCoord(transients.ra, transients.dec, frame='icrs', unit='deg')
        detection_skycoord = pixel_to_skycoord(detection.X_IMAGE, detection.Y_IMAGE, difference_wcs)
        transients_to_detection = truth_matching.skymatch_and_join(transients, detection, transients_skycoord, detection_skycoord, match_radius)
        detection_to_transients = truth_matching.skymatch_and_join(detection, transients, detection_skycoord, transients_skycoord, match_radius)
        
        transients_to_detection.to_csv(transients_to_detection_path, index=False)
        detection_to_transients.to_csv(detection_to_transients_path, index=False)
        return transients_to_detection, detection_to_transients
        
    def run(self):
        os.makedirs(self.base_output_dir, exist_ok=True)
        for i, row in self.data_records.iterrows():
            science_id = {
            'band': row['science_band'],
            'pointing': row['science_pointing'],
            'sca': row['science_sca']
            }
        
            template_id = {
            'band': row['template_band'],
            'pointing': row['template_pointing'],
            'sca': row['template_sca']
            }

            difference_id = self.__class__.get_difference_id(science_id, template_id)

            science_image_path = self.INPUT_IMAGE_PATTERN.format(**science_id)
            template_image_path = self.INPUT_IMAGE_PATTERN.format(**template_id)
            sub_output_dir = self.SUB_OUTPUT_DIR_PATTERN.format(**difference_id)
            full_output_dir = os.path.join(self.base_output_dir, sub_output_dir)
            os.makedirs(full_output_dir, exist_ok=True)

            # step 1: subtraction
            subtract = subtraction.Pipeline(science_band=science_id['band'],
                                            science_pointing=science_id['pointing'],
                                            science_sca=science_id['sca'],
                                            template_band=template_id['band'],
                                            template_pointing=template_id['pointing'],
                                            template_sca=template_id['sca'],
                                            out_dir=full_output_dir)
            subtract.run()
            
            # step 2: detection
            difference_image_path = os.path.join(full_output_dir, 'decorr_diff.fits')
            difference_detection_path = os.path.join(full_output_dir, 'detection.cat')
            source_detection.detect(difference_image_path, difference_detection_path,
                                    source_extractor_executable=self.SOURCE_EXTRACTOR_EXECUTABLE,
                                    detection_config=self.DETECTION_CONFIG
                                   )
            
            # step 3: truth retrival
            science_truth_path = self.INPUT_TRUTH_PATTERN.format(**science_id)
            template_truth_path = self.INPUT_TRUTH_PATTERN.format(**template_id)
            difference_truth_path = os.path.join(full_output_dir, 'truth.fits')

            truth = self.__class__.retrieve_truth(science_image_path, template_image_path,
                                   science_truth_path, template_truth_path,
                                   difference_truth_path)
            
            # step 4: truth matching
            transients_to_detection_path = os.path.join(full_output_dir, 'transients_to_detection.csv')
            detection_to_transients_path = os.path.join(full_output_dir, 'detection_to_transients.csv')
            transients_to_detection, detection_to_transients = self.__class__.match_transients(
                truth, difference_image_path, difference_detection_path, self.MATCH_RADIUS,
                transients_to_detection_path, detection_to_transients_path)
         
def main():
    parser = argparse.ArgumentParser( 'detection pipeline' )
    parser.add_argument( '-d', '--data_records', type=str, required=True, help="Input data records." )
    parser.add_argument( '-o', '--output_dir', type=str, default='../output', help="Output directory." )
    args = parser.parse_args()

    detection = Detection(args.data_records, args.output_dir)
    detection.run()

if __name__ == "__main__":
    main()
