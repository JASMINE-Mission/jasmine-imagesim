import numpy as np
import json
from jis.photonsim.wfe import wfe_model_z, calc_wfe, calc_dummy_wfe, calc_wfe_fringe37

def run_wfe(control_params, telescope):
    """ Making wfe. 
    
    Args:
        control_params: control parameters
        telescope: telescope object
    Returns:
        wavefront error

    """
    if control_params.effect.wfe == 'dummy':
        print("WFE simulation is skipped (making dummy).")
        wfe = calc_dummy_wfe(telescope.epd, 'dummy.json')
    elif control_params.effect.wfe == 'random':
        print("calculate random WFE...")
        wp  = control_params.wfe_control
        wfe_amplitudes = wfe_model_z(
            np.random, wp['zernike_nmax'], wp['reference_wl'],
            wp['zernike_odd'], wp['zernike_even'])

        # Saving amplitude data...
        with open(filenames["wfejson"], mode='w') as f:
            json.dump(wfe_amplitudes, f, indent=2)

        # Making wfe map...
        wfe = calc_wfe(telescope.epd, filenames["wfejson"])
    elif control_params.effect.wfe == 'fringe37':
        print("calculate WFE with fringe37 params...")
        wp = control_params.wfe_control

        # Making position array (in deg).
        positions = np.array([table_starplate['x pixel']-1.+detector.offset_x_mm/detector.pixsize/1.e-3,\
                              table_starplate['y pixel']-1.+detector.offset_y_mm/detector.pixsize/1.e-3]).T\
                    *detpix_scale/3600.
        ## detector.offset_[x|y]_mm is the position of (0, 0) on the telescope focal plane in mm.
        ## detector.pixsize is in um.

        # Making wfe map...
        wfe = calc_wfe_fringe37(telescope.epd, wp['fringe37_filename'],\
                                wp['reference_wl'], positions)
    else:
        print("WFE-calc method '{}' is not supported.".format(\
               control_params.effect.wfe))
        exit(-1)

    return wfe
