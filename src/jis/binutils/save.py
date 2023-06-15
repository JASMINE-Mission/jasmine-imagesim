import astropy.io.fits as pf
import numpy as np
import h5py
from jis.binutils.scales import get_pixelscales


def save_outputs(filenames, output_format, control_params, telescope, detector, wfe, psf, pixcube_global, tplate, uniform_flat_interpix, uniform_flat_intrapix, acex, acey, overwrite):
    """Saving the outputs."""
    # interpix flat
    if control_params.effect.flat_interpix is True:
        pf.writeto(filenames['interpix'],
                   detector.flat.interpix, overwrite=overwrite)
    else:
        pf.writeto(filenames['interpix'],
                   uniform_flat_interpix, overwrite=overwrite)

    # intrapix flat
    if control_params.effect.flat_interpix is True:
        pf.writeto(filenames['intrapix'],
                   detector.flat.intrapix, overwrite=overwrite)
    else:
        pf.writeto(filenames['intrapix'],
                   uniform_flat_intrapix, overwrite=overwrite)

    # psf
    hdu = pf.PrimaryHDU(psf)
    detpix_scale, fp_cellsize_rad, fp_scale, psfscale =
      get_pixelscales(control_params, telescope, detector)
    hdu.header['CDELT1'] = fp_scale 
    hdu.header['CDELT2'] = fp_scale
    hdu.header['CUNIT1'] = 'arcsec'
    hdu.header['CUNIT2'] = 'arcsec'
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filenames['psf'], overwrite=overwrite)

    # final image
    if output_format == 'hdfcube':
        with h5py.File(filenames['images'][0], 'w') as f:
            f.create_group('header')
            f.create_group('data')
            f.create_dataset('header/tplate', data=tplate)
            f.create_dataset('header/unit', data='e-/pix/plate')
            f.create_dataset('data/pixcube', data=pixcube_global)
    else:
        pixcube_global = np.swapaxes(pixcube_global, 0, 2)
        if output_format == 'platefits':
            for i in range(0, control_params.nplate):
                pf.writeto(filenames['images'][i], pixcube_global[i].astype(
                    'int32'), overwrite=overwrite)
        elif output_format == 'fitscube':
            pf.writeto(filenames['images'][0], pixcube_global.astype(
                'int32'), overwrite=overwrite)

    # wfe
    hdu = pf.PrimaryHDU(wfe)
    hdu.header['WFE-FILE'] = filenames['wfejson']
    hdu.header['WFE-EPD'] = telescope.epd
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filenames['wfe'], overwrite=overwrite)

    # aperture
    hdu = pf.PrimaryHDU(telescope.aperture)
    hdu.header['APTFILE'] = filenames['teljson']
    hdu.header['EPD'] = telescope.epd
    hdu.header['COBS'] = telescope.cobs
    hdu.header['STYPE'] = telescope.spider.type
    hdu.header['STEL'] = telescope.total_area  # total area in m^2
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filenames['aperture'], overwrite=overwrite)

    # ACE (x direction)
    hdu = pf.PrimaryHDU(acex)
    hdu.header['ACE-FILE'] = filenames['acejson']
    hdu.header['ACE-TOTT'] = control_params.ace_control['tace']
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filenames['acex'], overwrite=overwrite)

    # ACE (y direction)
    hdu = pf.PrimaryHDU(acey)
    hdu.header['ACE-FILE'] = filenames['acejson']
    hdu.header['ACE-TOTT'] = control_params.ace_control['tace']
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filenames['acey'], overwrite=overwrite)
