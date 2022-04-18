import json
import astropy.io.ascii as asc
from jis.photonsim.extract_json import Detector, ControlParams, Telescope


def load_parameters(filenames):
    """ Loading parameters. 

    Args:
        filenames: filename list

    Returns:
        table_starplate, detector, control_params, telescope, ace_params

    """
    table_starplate = asc.read(filenames["starplate"])
    detector        = Detector.from_json(filenames["detjson"])
    control_params  = ControlParams.from_json(filenames["ctljson"])
    telescope       = Telescope.from_json(filenames["teljson"])
    with open(filenames["acejson"], "r") as f:
        ace_params = json.load(f)
    f.close()
    return table_starplate, detector, control_params, telescope, ace_params
