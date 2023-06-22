

import inspect


class FunctionInspector:
    _FUNC_KWARGS = {}

    def __init__(self, func, ignore=[]):
        """
        Parameters:
        -----------
        func
            type: Any

        ignore
            type: list

        Returns:
        --------
        None
        """
        self._func = func
        self._ignore = ignore

    def _extract_func_params(self, func):
        return list(inspect.signature(func).parameters.keys())

    @property
    def FUNC_PARAMS(self):
        return self._extract_func_params(self._func)

    def forward(self, key, val):
        if (key in self.FUNC_PARAMS) and (not key in self._ignore):
            self._FUNC_KWARGS[key] = val

    def __call__(self, kwargs) -> dict:
        """
        Parameters:
        -----------
        kwargs
            type: dict

        Returns:
        --------
        func_kwargs
            type: dict
        """

        for key, val in kwargs.items():
            self.forward(key, val)

        return self._FUNC_KWARGS
