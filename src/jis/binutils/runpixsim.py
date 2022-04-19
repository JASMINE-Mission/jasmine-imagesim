import numpy as np
import json

def init_pix(control_params,detector, acex,acey, detpix_scale, driftsw):
    """ Preparation for making image, Setting and plotting full trajectory.

    Args:
        control_params: control parameters
        detector: detector object
        acex: ACE in x-axis
        acey: ACE in y-axis
        detpix_scale: scale of detector pixels 
        driftsw: bool if True drift included

    Returns:
        theta_full, pixdim, Npixcube

    """

    if driftsw:
        from jis.photonsim.extract_json import Drift
        dft=Drift.from_json(filenames["dftjson"])
        dft.compute_drift(control_params.ace_control['dtace'],len(acex))
        
    ## Full data of the displacement in detpix.
    ## (ace[x|y] scaled and converted to detpix)
    acex_std = control_params.ace_control.get('acex_std')
    acey_std = control_params.ace_control.get('acey_std')
    if driftsw:
        theta_full = np.array([acex*acex_std/detpix_scale+dft.drift_theta[0,:], acey*acey_std/detpix_scale+dft.drift_theta[1,:]])
        plt.plot(acex*acex_std/detpix_scale+dft.drift_theta[0,:], acey*acey_std/detpix_scale+dft.drift_theta[1,:], ".")
        plt.savefig("theta.png")
    else:
        theta_full = np.array([acex*acex_std/detpix_scale, acey*acey_std/detpix_scale])

    Npixcube = int((np.max(np.abs(theta_full))+detector.nmargin)*2)
    pixdim   = [Npixcube, Npixcube] # adaptive pixel dimension in the aperture.
    return theta_full, pixdim, Npixcube

