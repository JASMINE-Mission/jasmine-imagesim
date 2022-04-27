import os


def set_output_from_args(args):
    """set output from docopt object.

    Args:
        args: docopt(__doc__)

    Returns:
        output_format, overwrite
    """
    output_format = args['--format']
    if output_format not in ['platefits', 'fitscube', 'hdfcube']:
        print("format must be 'platefits', 'fitscube' or 'hdfcube'.")
        exit(-1)

    overwrite = False
    if args['--overwrite']:
        overwrite = True

    return output_format, overwrite


def check_output_directory(filenames, dirname_output, overwrite):
    """Checking the output directory.

    Args:
       filenames: filename to be checked
       dirname_output: output direction
       overwrite: overwrite status

    """
    if not os.path.exists(dirname_output):
        os.makedirs(dirname_output)
    else:
        if overwrite is not True:
            for filename in filenames['output']:
                if os.path.exists(filename):
                    print("\"{}\" exists.".format(filename))
                    print('Please set --overwrite option to overwrite it.')
                    exit(-1)


def set_filenames_from_args(args):
    """
    Args:
        args: docopt(__doc__)

    Returns:
        filename dictionary, output directory

    """
    # Getting the parameters from command line. ####################
    dirname_params = ''
    if args['--pd']:
        dirname_params = args['--pd']

    filenames = {}
    filenames['starplate'] = os.path.join(dirname_params, args['--starplate'])
    if args['--var']:
        filenames['varjson'] = os.path.join(dirname_params, args['--var'])
    filenames['detjson'] = os.path.join(dirname_params, args['--det'])
    filenames['teljson'] = os.path.join(dirname_params, args['--tel'])
    filenames['acejson'] = os.path.join(dirname_params, args['--ace'])
    filenames['ctljson'] = os.path.join(dirname_params, args['--ctl'])
    if args['--dft']:
        filenames['dftjson'] = os.path.join(dirname_params, args['--dft'])

    # Setting output filenames. ####################################
    dirname_output = '.'
    if args['--od']:
        dirname_output = args['--od']
    filenames['interpix'] = os.path.join(dirname_output, 'interpix.fits')
    filenames['intrapix'] = os.path.join(dirname_output, 'intrapix.fits')
    filenames['wfejson'] = os.path.join(dirname_output, 'wfe.json')
    filenames['wfe'] = os.path.join(dirname_output, 'wfe.fits')
    filenames['aperture'] = os.path.join(dirname_output, 'aperture.fits')
    filenames['psf'] = os.path.join(dirname_output, 'psf.fits')
    filenames['acex'] = os.path.join(dirname_output, 'aceX.fits')
    filenames['acey'] = os.path.join(dirname_output, 'aceY.fits')

    return filenames, dirname_output


def set_filenames_output(args, filenames, control_params, dirname_output):
    """set output.

    Args:
        args: docopt(__doc__)
        filenames: filenames list
        control_params: control parameter
        dirname_output: output directory name

    Returns:
        filenames, output_format, overwrite
    """
    output_format, overwrite = set_output_from_args(args)
    if output_format == 'platefits':
        filenames['images'] = []
        for i in range(0, control_params.nplate):
            filenames['images'].append(os.path.join(
                dirname_output, 'image{:02d}.fits'.format(i)))
    elif output_format == 'fitscube':
        filenames['images'] = [os.path.join(dirname_output, 'image.fits')]
    elif output_format == 'hdfcube':
        filenames['images'] = [os.path.join(dirname_output, 'image.h5')]

    filenames['output'] = [filenames['interpix'], filenames['intrapix'],
                           filenames['wfejson'], filenames['wfe'],
                           filenames['aperture'], filenames['psf'],
                           filenames['acex'], filenames['acey']]
    filenames['output'] = filenames['output'] + filenames['images']

    check_output_directory(filenames, dirname_output, overwrite)
    return filenames, output_format, overwrite
