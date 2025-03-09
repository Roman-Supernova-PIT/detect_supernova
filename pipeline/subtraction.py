import io
import shutil
import os
import pathlib
import argparse
from dataclasses import dataclass
from astropy.io import fits
import numpy as np
import cupy as cp
from sfft.SpaceSFFTCupyFlow import SpaceSFFT_CupyFlow
from roman_imsim.utils import roman_utils


INPUT_IMAGE_PATTERN = ("RomanTDS/images/simple_model/{band}/{pointing}/Roman_TDS_simple_model_{band}_{pointing}_{sca}.fits.gz")
INPUT_TRUTH_PATTERN = ("RomanTDS/truth/{band}/{pointing}/Roman_TDS_index_{band}_{pointing}_{sca}.txt")
SIMS_DIR = pathlib.Path( os.getenv( 'SIMS_DIR', None ) )
GALSIM_CONFIG = pathlib.Path( os.getenv("SN_INFO_DIR" ) ) / "tds.yaml"

IMAGE_WIDTH =4088
IMAGE_HEIGHT = 4088


def gz_and_ext(in_path,out_path):
    # Modified from https://github.com/Roman-Supernova-PIT/phrosty/blob/main/phrosty/imagesubtraction.py#L77
    
    """
    Utility function that unzips the original file and turns it into a single-extension FITS file.
    """

    bio = io.BytesIO()
    with gzip.open(in_path,'rb') as f_in:
        shutil.copyfileobj(f_in, bio)
    bio.seek(0)

    with fits.open(bio) as hdu:
        newhdu = fits.HDUList([fits.PrimaryHDU(data=hdu[1].data, header=hdu[0].header)])
        newhdu.writeto(out_path, overwrite=True)

    return out_path

def sky_subtract( inpath, skysubpath, detmaskpath, temp_dir=pathlib.Path("/tmp"), force=False ):
    # Modified from https://github.com/Roman-Supernova-PIT/phrosty/blob/main/phrosty/imagesubtraction.py#L100
    
    """Subtracts background, found with Source Extractor.

    Parameters
    ----------
      inpath: Path
        Original FITS image

      skysubpath: Path
        Sky-subtracted FITS image

      detmaskpath: Path
        Detection Mask FITS Image.  (Will be uint8, I think.)

      temp_dir: Path
        Already-existing directory where we can write a temporary file.
        (If the image is .gz compressed, source-extractor can't handle
        that, so we have to write a decompressed version.)

      force: bool, default False
        If False, and outpath already exists, do nothing.  If True,
        clobber the existing file and recalculate it.

    Returns
    -------
      skyrms: float
        Median of the skyrms image calculated by source-extractor

    """

    if ( not force ) and ( skysubpath.is_file() ) and ( detmaskpath.is_file() ):
        with fits.open( skysubpath ) as hdul:
            skyrms = hdul[0].header['SKYRMS']
        return skyrms

    if inpath.name[-3:] == '.gz':
        decompressed_path = temp_dir / inpath.name[:-3]
        gz_and_ext( inpath, decompressed_path )
    else:
        decompressed_path = inpath


    ( SKYDIP, SKYPEAK, PixA_skysub,
      PixA_sky, PixA_skyrms ) = SEx_SkySubtract.SSS(FITS_obj=decompressed_path,
                                                    FITS_skysub=skysubpath,
                                                    FITS_detmask=detmaskpath,
                                                    FITS_sky=None, FITS_skyrms=None,
                                                    ESATUR_KEY='ESATUR',
                                                    BACK_SIZE=64, BACK_FILTERSIZE=3,
                                                    DETECT_THRESH=1.5, DETECT_MINAREA=5,
                                                    DETECT_MAXAREA=0,
                                                    VERBOSE_LEVEL=2, MDIR=None)

    return np.median( PixA_skyrms )

def get_imsim_psf(x, y, pointing, sca, size=201, config_yaml_file=None,
                  psf_path=None, **kwargs):
    # Modified from https://github.com/Roman-Supernova-PIT/phrosty/blob/main/phrosty/imagesubtraction.py#L239C1-L241C1
    # It seems that band information is not needed

    config_path = config_yaml_file
    config = roman_utils(config_path, pointing, sca)
    psf = config.getPSF_Image(size, x, y, **kwargs)
    psf.write(str(psf_path))

@dataclass
class ImageInfo:
    data_id: dict
    
    def __post_init__(self):
        self.image_path = SIMS_DIR / pathlib.Path(INPUT_IMAGE_PATTERN.format(**self.data_id))
        self.cx = IMAGE_WIDTH // 2
        self.cy = IMAGE_HEIGHT // 2

        self.skysub_path = None
        self.detect_mask_path = None
        self.skyrms = None
        self.psf_path = None

        
class Pipeline:

    def __init__(self, science_band, science_pointing, science_sca,
                 template_band, template_pointing, template_sca,
                 galsim_config_file=GALSIM_CONFIG, temp_dir="/phrosty_temp", out_dir="./output"):

        self.science_info = ImageInfo({'band': science_band, 'pointing': science_pointing, 'sca': science_sca})
        self.template_info = ImageInfo({'band': template_band, 'pointing': template_pointing, 'sca': template_sca})
        
        self.galsim_config_file = galsim_config_file
        self.temp_dir = pathlib.Path(temp_dir)
        self.out_dir = pathlib.Path(out_dir)
        os.makedirs(self.out_dir, exist_ok=True)

    def run_sky_subtract(self, image_info):
        image_name = image_info.image_path.name
        skysub_path = self.temp_dir / f"skysub_{image_name}"
        detmask_path = self.temp_dir / f"detmask_{image_name}"
        skyrms = sky_subtract(image_info.image_path, skysub_path, detmask_path, temp_dir=self.temp_dir, force=False )
        
        image_info.skysub_path = skysub_path
        image_info.detmask_path = detmask_path
        image_info.skyrms = skyrms

    def run_get_imsim_psf(self, image_info):
        image_name = image_info.image_path.name
        psf_path = self.temp_dir / f"psf_{image_name}"
        
        get_imsim_psf(x=image_info.cx, y=image_info.cy,
                      pointing=image_info.data_id['pointing'], sca=image_info.data_id['sca'],
                      size=201, psf_path=psf_path, config_yaml_file=self.galsim_config_file, include_photonOps=True)

        image_info.psf_path = psf_path

    def align_and_pre_convolve(self):

        with fits.open( self.science_info.skysub_path ) as hdul:
            hdr_sci = hdul[0].header
            data_sci = cp.array( np.ascontiguousarray(hdul[0].data.T), dtype=cp.float64 )
            
        with fits.open( self.template_info.skysub_path ) as hdul:
            hdr_templ = hdul[0].header
            data_templ = cp.array( np.ascontiguousarray(hdul[0].data.T), dtype=cp.float64 )

        with fits.open( self.science_info.psf_path ) as hdul:
            sci_psf = cp.array( np.ascontiguousarray( hdul[0].data.T ), dtype=cp.float64 )

        with fits.open( self.template_info.psf_path ) as hdul:
            templ_psf = cp.array( np.ascontiguousarray( hdul[0].data.T ), dtype=cp.float64 )

        with fits.open( self.science_info.detmask_path ) as hdul:
            sci_detmask = cp.array( np.ascontiguousarray( hdul[0].data.T ) )

        with fits.open( self.template_info.detmask_path ) as hdul:
            templ_detmask = cp.array( np.ascontiguousarray( hdul[0].data.T ) )

        sfftifier = SpaceSFFT_CupyFlow(
            hdr_sci, hdr_templ,
            self.science_info.skyrms,
            self.template_info.skyrms,
            data_sci, data_templ,
            sci_detmask, templ_detmask,
            sci_psf, templ_psf
        )
        sfftifier.resampling_image_mask_psf()
        sfftifier.cross_convolution()

        return sfftifier

    def run_decorrelation(self, sfftifier):

        decorr_diff = sfftifier.apply_decorrelation(sfftifier.PixA_DIFF_GPU)
        decorr_zptimg = sfftifier.apply_decorrelation(sfftifier.PixA_Ctarget_GPU)
        decorr_psf = sfftifier.apply_decorrelation(sfftifier.PSF_target_GPU)

        decoor_diff_path = self.out_dir / "decorr_diff.fits"
        decoor_zptimg_path = self.out_dir / "decorr_zptimg.fits"
        decorr_psf_path = self.out_dir / "decorr_psf.fits"

        fits.writeto( decoor_diff_path, cp.asnumpy(decorr_diff).T, header=sfftifier.hdr_target, overwrite=True )
        fits.writeto( decoor_zptimg_path, cp.asnumpy(decorr_zptimg).T, header=sfftifier.hdr_target, overwrite=True )
        fits.writeto( decorr_psf_path, cp.asnumpy(decorr_psf).T, header=None, overwrite=True )

    def run(self):
        # sky subtraction
        self.run_sky_subtract(self.science_info)
        self.run_sky_subtract(self.template_info)

        # get psf
        self.run_get_imsim_psf(self.science_info)
        self.run_get_imsim_psf(self.template_info)
        
        # prec-onvolution
        sfftifier = self.align_and_pre_convolve()

        # subtraction
        sfftifier.sfft_subtraction()

        # decorrelation
        sfftifier.find_decorrelation()
        self.run_decorrelation(sfftifier)    
        
def main():
    parser = argparse.ArgumentParser( 'phrosty pipeline' )
    parser.add_argument('--science_band', type=str, required=True, help="Science band" )
    parser.add_argument('--science_pointing', type=int, required=True, help="Science pointing")
    parser.add_argument('--science_sca', type=int, required=True, help="Science sca")
    parser.add_argument('--template_band', type=str, required=True, help="Template band" )
    parser.add_argument('--template_pointing', type=int, required=True, help="Template pointing")
    parser.add_argument('--template_sca', type=int, required=True, help="Template sca")
    # parser.add_argument( '-p', '--nprocs', type=int, default=1, help="Number of process for multiprocessing steps (e.g. skysub)" )
    # parser.add_argument( '-w', '--nwrite', type=int, default=5, help="Number of parallel FITS writing processes" )
    # parser.add_argument( '-v', '--verbose', action='store_true', default=False, help="Show debug log info" )
    parser.add_argument( '--out-dir', default="/out_dir", help="Output dir, default /out_dir" )
    # parser.add_argument( '--ltcv-dir', default="/lc_out_dir", help="Output dir for lightcurves, default /lc_out_dir" )
    parser.add_argument( '--temp-dir', default="/phrosty_temp", help="Temporary working dir, default /phrosty_temp" )


    args = parser.parse_args()

    galsim_config = pathlib.Path( os.getenv("SN_INFO_DIR" ) ) / "tds.yaml"

    pipeline = Pipeline(args.science_band, args.science_pointing, args.science_sca,
                        args.template_band, args.template_pointing, args.template_sca,
                        galsim_config_file=galsim_config,
                        temp_dir=args.temp_dir, out_dir=args.out_dir)
    """
    pipeline = Pipeline( args.ra, args.dec, args.band, science_images, template_images,
                         nprocs=args.nprocs, nwrite=args.nwrite,
                         temp_dir=args.temp_dir, out_dir=args.out_dir, ltcv_dir=args.ltcv_dir,
                         galsim_config_file=galsim_config, force_sky_subtract=args.force_sky_subtract,
                         nuke_temp_dir=False, verbose=args.verbose )
    """
    pipeline.run()

# ======================================================================
if __name__ == "__main__":
    main()