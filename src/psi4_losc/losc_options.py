"""
Options to configure LOSC calculations within psi4_losc package.
"""


class Options():
    """Options class that controls the LOSC calculations in `psi4_losc` module.

    The options are specified by three labels:

        - `module`: the module of the options.
        - `key`: the key of the option in `module`.
        - `value`: the value of the option in `module`.

    All `module` and `key` values are str and case-insensitive.

    - **The curvature module**: `module = "curvature"`

      Below lists all the valid `key`-`value` options for the curvature module.

        - version : {2, 1}, default=2.
            The version of LOSC curvature.
        - v1_parameter_tau: float, default=1.2378
            The parameter :math:`\\tau` in LOSC curvature version 1.
            Require `version=1` to use this setting.
        - v2_parameter_tau: float, default=1.2378
              The parameter :math:`\\tau` in LOSC curvature version 2.
              Require `version=2` to use this setting.
        - v2_parameter_zeta: float, default=8.0
              The parameter :math:`\\zeta` in LOSC curvature version 2.
              Require `version=2` to use this setting.
        - df_molecular_fragment_size: int, default=2
              The size in the number of atoms to split the module. This is used
              to achieve the block-wise construction of two-electron integral of
              curvature with density fitting.

    - **The localizer module**: `module = "localizer"`

      Below lists all the valid `key`-`value` options for the localization module.

        - version : {2}, default=2.
            The version of LOSC localizer.
        - max_iter: int, default=1000
            The maximum number of iterations in localization.
        - convergence: float, default=1e-10
            The convergence tolerance for the localization.
        - random_permutation: bool, default=True
            Use the random permutation for Jacob-Sweep algorithm in the
            localization or not.
        - v2_parameter_gamma: float, default=0.707
            The parameter :math:`\\gamma` in LOSC localizer version 2.
            Require `version=2` to use this setting.
        - v2_parameter_c: float, default=1000
            The parameter :math:`C` in LOSC localizer version 2.
            Require `version=2` to use this setting.
    """

    def __init__(self):
        self._options = {
            'curvature': {
                'version': 2,
                'v1_parameter_tau': 1.2378,
                'v2_parameter_tau': 1.2378,
                'v2_parameter_zeta': 8.0,
                'df_molecular_fragment_size': 2,
            },
            'localizer': {
                'version': 2,
                'max_iter': 1000,
                'convergence': 1e-10,
                'random_permutation': True,
                'v2_parameter_gamma': 0.707,
                'v2_parameter_c': 1000,
            }
        }

    def set_param(self, module, key, value):
        """Set an option for the LOSC calculation.

        Parameters
        ----------
        module : str
            The module name of the option.
        key : str
            The key name of the option in option module `module`.
        value :
            The value name of the option in option module `module`.

        See also
        --------
        Options : all valid options.
        """
        module, key = map(str.lower, (module, key))
        if module not in self._options:
            raise Exception(
                f'Wrong module name of LOSC options: module="{module}"')
        if key not in self._options[module]:
            raise Exception(
                f'Wrong key to set options of LOSC: module="{module}" key={key}')
        self._options[module][key] = value

    def set_params(self, options):
        """Set one or more options for the LOSC calculation.

        Parameters
        ----------
        options : dict {module: {key: value}}
            The options dictionary.

        See also
        --------
        Options : all valid options.
        """

        for module, opts in options.items():
            module = module.lower()
            if module not in self._options:
                raise Exception(
                    f'Wrong module name of LOSC options: module="{module}"')
            for key, val in opts.items():
                key = key.lower()
                if key not in self._options[module]:
                    raise Exception(
                        f'Wrong key to set options of LOSC: module="{module}" key={key}')
                self._options[module][key] = val

    def get_param(self, module, key):
        """Get the value of an option.

        Parameters
        ----------
        module : str
            The module name of the option.
        key : str
            The key name of the option in option module `module`.

        Returns
        -------
        value:
            The value of the option.

        See also
        --------
        Options : all valid options and the return type of the value.
        """
        module, key = map(str.lower, (module, key))
        try:
            rst = self._options[module][key]
        except Exception:
            raise Exception(
                f'Unknown module or key to get the LOSC option: module="{module}" key={key}')
        return rst

    def __str__(self):
        """Print all LOSC options to psi4 output file."""
        rst = []
        rst.append('')
        rst.append('==> LOSC Settings <==')
        for module, opts in self._options.items():
            rst.append(f'--> module={module}')
            for key, val in opts.items():
                rst.append(f'    {key}: {val}')
            rst.append('')
        return '\n'.join(rst)
