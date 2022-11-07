import numpy as np

def load_filenames_for_test():
    import pkg_resources
    import os
    dirname_params = pkg_resources.resource_filename(
        'jis', 'data/params')
    filenames={}
    filenames['starplate'] = os.path.join(dirname_params, "star_plate.csv")
    filenames['detjson'] = os.path.join(dirname_params, "det.json")
    filenames['teljson'] = os.path.join(dirname_params, "tel.json")
    filenames['acejson'] = os.path.join(dirname_params, "ace_001.json")
    filenames['ctljson'] = os.path.join(dirname_params, "ctl.json")
    return filenames

def pixel_scale_radian(detector,telescope):
    return detector.pixsize*1.e-6/(telescope.efl*1.e-3)

def test_pixel_scale_radian(detector,telescope):
    filenames = load_filenames_for_test()
    table_starplate, detector, control_params, telescope, ace_params = load_parameters(
        filenames)
    pixel_scale_radian(detector,telescope)

def set_wcs(detector,telescope):
    from astropy import wcs 
    w = wcs.WCS(naxis=2)
    w.wcs._naxis = [detector.npix, detector.npix] #the number of the one side pixels 
    w.wcs.crpix = [detector.npix/2.0, detector.npix/2.0] #reference pixel coordinate
    w.wcs.cdelt = [pixel_scale_radian(detector,telescope), pixel_scale_radian(detector,telescope)]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    pixel_scale_radian = 
  
def test_set_wcs():
    from jis.binutils.setcontrol import load_parameters
    import os
    
    filenames = load_filenames_for_test()
    table_starplate, detector, control_params, telescope, ace_params = load_parameters(
        filenames)
    print()
    set_wcs(detector,telescope)


if __name__ == "__main__":
    
    test_set_wcs()