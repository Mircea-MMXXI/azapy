from .CVaRAnalyzer import CVaRAnalyzer
from .Port_InvVol import Port_InvVol


class Port_CVaR(Port_InvVol):
    """
    Backtesting mCVaR optimal portfolio strategies, periodically rebalanced.

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

    def set_model(self, alpha=[0.975], coef=None, rtype='Sharpe',
                  mu=None, mu0=0, aversion=None, ww0=None, 
                  hlength=3.25, method='ecos'):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        `alpha` : list, optional
            List of alpha confidence levels. The default is `[0.975]`.
        `coef` : list, optional
            List of positive mixture coefficients. Note that `len(coef)`
            must be equal to `len(alpha)`. A `None` value assumes an
            equal weighted risk mixture.
            The vector of coefficients will be normalized to unit.
            The default is `None`.
        `rtype` : str, optional
            Optimization type. Possible values \n
                'Risk' : minimization of dispersion (risk) measure for  
                targeted rate of return. \n
                'Sharpe' : maximization of generalized Sharpe ratio.\n
                'Sharpe2' : minimization of the inverse generalized Sharpe 
                ratio.\n
                'MinRisk' : minimum dispersion (risk) portfolio.\n
                'InvNrisk' : optimal portfolio with the same dispersion (risk)
                value as a benchmark portfolio 
                (e.g. equal weighted portfolio).\n
                'RiskAverse' : optimal portfolio for a fixed value of 
                risk-aversion factor.
            The default is `'Sharpe'`.
        `mu` : float, optional
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype='Risk'`
            The default is `None`.
        `mu0` : float, optional
            Risk-free rate accessible to the investor.
            Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        `aversion` : float, optional
            The value of the risk-aversion factor.
            Must be positive. Relevant only if `rtype='RiskAvers'`.
            The default is `None`.
        `ww0` : list (also `numpy.array` or `pandas.Series`), optional
            Targeted portfolio weights. 
            Relevant only if `rype='InvNrisk'`.
            Its length must be equal to the number of
            symbols in rrate (mktdata). 
            All weights must be >= 0 with sum > 0.
            If it is a list or a `numpy.array` then the weights are assumed to
            by in order of `rrate.columns`. If it is a `pandas.Series` then 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (same symbols, not necessary in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.
        `hlength` : float, optional
            The length in year of the historical calibration period relative
            to 'Dfix'. A fractional number will be rounded to an integer number
            of months. The default is `3.25` years.
        `method` : str, optional
            Linear programming numerical method.
            Could be: 'ecos', 'highs-ds', 'highs-ipm', 'highs',
            'interior-point', 'glpk' and 'cvxopt'.
            The default is `'ecos'`.

         Returns
        -------
        `pandas.DataFrame`
            The portfolio time-series in the format 'date', 'pcolname'.
        """
        self._set_alpha(alpha, coef)
        self._set_rtype(rtype)
        self.mu = mu
        self.mu0 = mu0
        self.aversion = aversion
        self.ww0 = ww0
        self.hlength = hlength
        self._set_method(method)

        self._set_schedule()
        self._set_weights()
        self._port_calc()
        return self.port


    def _set_alpha(self, alpha, coef):
        # alpha
        self.alpha = alpha
        self.coef = coef


    def _set_rtype(self, rtype):
        rtype_values = ['Sharpe', 'Risk', 'MinRisk', 'InvNrisk', 'RiskAverse',
                        'Sharpe2']
        if not rtype in rtype_values:
            raise ValueError(f"rtype must be one of {rtype_values}")

        self.rtype = rtype


    def _set_method(self, method):
        self.method = method


    def _wwgen(self):
        return CVaRAnalyzer(self.alpha, self.coef, rtype=self.rtype,
                            method=self.method)


    def _ww_calc(self, data):
        return self._wwgen().getWeights(mu=self.mu, mu0=self.mu0,
                aversion=self.aversion, ww0 = self.ww0, rrate=data)
