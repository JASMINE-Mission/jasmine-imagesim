def load_filenames_for_test():
    """load params for test

    Returns:
        detector, control_params, telescope instances
    """
    import pkg_resources
    import os
    from jis.photonsim.extract_json import Detector, ControlParams, Telescope
    dirname_params = pkg_resources.resource_filename('jis', 'data/params')
    filenames = {}
    filenames['detjson'] = os.path.join(dirname_params, "det.json")
    filenames['teljson'] = os.path.join(dirname_params, "tel.json")
    filenames['ctljson'] = os.path.join(dirname_params, "ctl.json")

    detector = Detector.from_json(filenames['detjson'])
    control_params = ControlParams.from_json(filenames['ctljson'])
    telescope = Telescope.from_json(filenames['teljson'])

    return detector, control_params, telescope
