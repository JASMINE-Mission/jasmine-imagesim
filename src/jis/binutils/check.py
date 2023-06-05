def check_ace_length(Nts_per_plate, control_params, theta_full):
    """check the length of ace
    Args:
        Nts_per_plate: number of plates
        control_params: contorl parameters
        theta_full: full trajectory

    """

    if Nts_per_plate*control_params.nplate >= theta_full.shape[1]:
        print('Insufficient time length of ACE data.')
        print('Nts_per_plate*Nplate: {}'.format(Nts_per_plate*control_params.nplate))
        print('N_ace_data          : {}'.format(theta_full.shape[1]))
        raise ValueError('error!')
    else:
        print('ACE length check passed.')
