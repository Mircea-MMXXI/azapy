import numpy as np

from .CVaRAnalyzer import CVaRAnalyzer
from .Port_InvVol import Port_InvVol


class Port_CVaR(Port_InvVol):
    """
    Back testing the CVaR optimal portfolio periodically rebalanced.

    Methods:
        * set_model
        * get_port
        * get_nshares
        * get_weights
        * get_account
        * get_mktdata
        * port_view
        * port_view_all
        * port_drawdown
        * port_perf
        * port_annual_returns
        * port_monthly_returns
        * port_period_returns
    """
    def set_model(self, mu, alpha=[0.975], coef=None, rtype='Sharpe',
                  hlength=3.25, method='ecos'):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        mu : float
            Reference rate. Its meaning depends on the value of rtype. For
            rtype equal to: \n
                'Sharpe' : `mu` is the risk-free rate. \n
                'Ris' : `mu` is the targeted expected rate of returns. \n
                'MinRisk' and 'InvNrisk' : `mu` is ignored. \n
                'RiskAverse' : `mu` is the Lambda risk aversion coefficient.
        alpha : list, optional
            List of alpha confidence levels. The default is [0.975].
        coef : list, optional
            List of mixture coefficients values. Note that `len(coef)`
            must be equal to `len(alpha)`. A value of `None` assumes
            `coef = [1 / len(alpha)] * len(alpha)`.
       rtype : str, optional
            Optimization type. Possible values \n
                "Risk" : minimization of dispersion (risk) measure for a fixed 
                vale of expected rate of return. \n
                "Sharpe" : maximization of generalized Sharpe ratio.\n
                "Sharpe2" : minimization of the inverse generalized Sharpe 
                ratio.\n
                "MinRisk" : optimal portfolio with minimum dispersion (risk) 
                value.\n
                "InvNRisk" : optimal portfolio with the same dispersion (risk)
                value as equal weighted portfolio. \n
                "RiskAverse" : optimal portfolio for a fixed value of risk 
                aversion coefficient.
            The default is "Sharpe".
        hlength : float, optional
            The length in year of the historical calibration period relative
            to 'Dfix'. A fractional number will be rounded to an integer number
            of months. The default is 3.25 years.
        method : str, optional
            Linear programming numerical method.
            Could be: 'ecos', 'highs-ds', 'highs-ipm', 'highs',
            'interior-point', 'glpk' and 'cvxopt'.
            The default is 'ecos'.

        Returns
        -------
        pandas.DataFrame
            The portfolio time-series in the format 'date', 'pcolname'.
        """
        self._set_alpha(alpha, coef)
        self._set_rtype(rtype)
        self.mu = mu
        self.hlength = hlength
        self._set_method(method)

        self._set_schedule()
        self._set_weights()
        self._port_calc()
        return self.port

    def _set_alpha(self, alpha, coef):
        # alpha
        self.alpha = np.array(alpha)
        if np.any((self.alpha <= 0.) | (1. <= self.alpha)):
            raise ValueError("alpha must be in (0, 1)")

        # coef
        if coef is None:
            self.coef = np.ones(len(self.alpha))
        else:
            if len(coef) != len(self.alpha):
                raise ValueError("coef must have same length as alpha")
            self.coef = np.array(coef)

        if np.any(self.alpha < 0.):
            raise ValueError("coef must be >= 0")

        scoef = self.coef.sum()
        if scoef <= 0.:
            raise ValueError("at leas one coef must be > 0")

        self.coef = self.coef / scoef

    def _set_rtype(self, rtype):
        rtype_values = ['Sharpe', 'Risk', 'MinRisk', 'InvNrisk', 'RiskAverse',
                        'Sharpe2']
        if not rtype in rtype_values:
            raise ValueError(f"rtype must be one of {rtype_values}")

        self.rtype = rtype

    def _set_method(self, method):
        method_values = ['ecos', 'highs-ds', 'highs-ipm', 'highs',
                         'interior-point', 'glpk', 'cvxopt']
        if not method in method_values:
            raise ValueError(f"mehtod must be one of {method_values}")
        self.method = method

    def _wwgen(self):
        return CVaRAnalyzer(self.alpha, self.coef, rtype=self.rtype,
                            method=self.method)

    def _ww_calc(self, data):
        return self._wwgen().getWeights(mu=self.mu, rrate=data)
