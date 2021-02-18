"""
Options to configure LOSC calculations within psi4_losc package.
"""

from psi4_losc import utils

# ==> The LOSC curvature options <==
curvature = {
    'version': 2,
}


# ==> The LOSC localization options <==
localization = {
    'version': 2,
    'max_iter': 1000,
    'convergence': 1e-10,
    'random_permutation': True,
    'v2_parameter_gamma': 0.707,
    'v2_parameter_c': 1000,
}


def show_options():
    local_print = utils.init_local_print(float('inf'))
    local_print(0, '')
    local_print(0, '==> LOSC Curvature Settings <==')
    local_print(0, 'version: ', curvature['version'])
    local_print(0, '')

    local_print(0, '==> LOSC Localization Settings <==')
    local_print(0, 'version: ', localization['version'])
    local_print(0, 'max_iter: ', localization['max_iter'])
    local_print(0, 'random_permutation: ', localization['random_permutation'])
    if localization['version'] == 2:
        local_print(0, 'parameter gamma: ', localization['v2_parameter_gamma'])
        local_print(0, 'parameter C: ', localization['v2_parameter_c'])
    local_print(0, '')
