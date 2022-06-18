from .Port_CVaR import Port_CVaR
from .BTADAnalyzer import BTADAnalyzer

class Port_BTAD(Port_CVaR):
    """
    Backtesting mBTAD optimal portfolio strategies, periodically rebalanced.

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
    def set_model(self, alpha=[0.], coef=None, rtype='Sharpe', 
                  mu=None, mu0=0., aversion=None, ww0=None,
                  detrended=False, hlength=3.25, method='ecos'):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        `alpha` : list, optional
            List of thresholds. The default is `[0.]`.
        `coef` : list, optional
            List of positive mixture 
            coefficients. Must have the same size as `alpha`. 
            A `None` value assumes an equal weighted risk mixture.
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
        `detrended` : Boolean, optional
            Designates the rate type used in the delta-risk calculations:\n
                `True` : detrended rate of return, i.e. r - E(r), \n
                `False` : standard rate of return. 
            The default is `False`.
        `hlength` : float, optional
            The length in year of the historical calibration period relative
            to 'Dfix'. A fractional number will be rounded to an integer number
            of months. The default is `3.25` years.
        `method` : str, optional
            Linear programming numerical method.
            Could be: 'ecos', 'highs-ds', 'highs-ipm', 'highs',
            'interior-point', 'glpk' and 'cvxopt'.
            The defualt is `'ecos'`.

        Returns
        -------
        pandas.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        self.detrended = detrended
        return super().set_model(alpha=alpha, coef=coef, rtype=rtype, 
                                 mu=mu, mu0=mu0, aversion=aversion, ww0=ww0,
                                 hlength=hlength, method=method)


    def _wwgen(self):
        return BTADAnalyzer(self.alpha, self.coef, rtype=self.rtype,
                            detrended=self.detrended, method=self.method)
