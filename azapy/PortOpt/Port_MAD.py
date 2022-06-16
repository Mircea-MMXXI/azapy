from .Port_CVaR import Port_CVaR
from .MADAnalyzer import MADAnalyzer


class Port_MAD(Port_CVaR):
    """
    Back testing the MAD optimal portfolio periodically rebalanced.

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
    def set_model(self, coef=[1.], rtype='Sharpe', mu=None, mu0=0,
                  aversion=None, ww0=None, hlength=3.25, method='ecos'):
        """
        Set model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        `coef` : list, optional
            Positive non-increasing list of mixture coefficients. 
            The default is [1.].
        `rtype` : : str, optional
            Optimization type. Possible values \n
                'Risk' : minimization of dispersion (risk) measure for a fixed 
                vale of expected rate of return. \n
                'Sharpe' : maximization of generalized Sharpe ratio.\n
                'Sharpe2' : minimization of the inverse generalized Sharpe 
                ratio.\n
                'MinRisk' : optimal portfolio with minimum dispersion (risk) 
                value.\n
                'InvNRisk' : optimal portfolio with the same dispersion (risk)
                as the targeted portfolio
                (e.g. equal weighted portfolio). \n
                'RiskAverse' : optimal portfolio for a fixed value of risk 
                aversion coefficient.
            The default is 'Sharpe'.
        `mu` : float, optional
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype='Risk'`
            The default is `None`.
        `mu0` : float, optional
            Risk-free rate accessible to the investor.
            Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        `aversion` : float, optional
            The value of the risk-aversion coefficient.
            Must be positive. Relevant only if `rtype='RiskAvers'`.
            The default is `None`.
        `ww0` : list (also np.array or pandas.Series), optional
            Targeted portfolio weights. 
            Relevant only if rype='InvNrisk'.
            Its length must be equal to the number of
            symbols in rrate (mktdata). 
            All weights must be >= 0 with sum > 0.
            If it is a list or a numpy.array then the weights are assumed to
            by in order of rrate.columns. If it is a pandas.Series the index
            should be compatible with the rrate.columns or mktdata symbols
            (same symbols, not necessary in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.
        `hlength` : float, optional
            The length in year of the historical calibration period relative
            to 'Dfix'. A fractional number will be rounded to an integer number
            of months. The default is 3.25 years.
        `method` : str, optional
            Linear programming numerical method.
            Could be: 'ecos', 'highs-ds', 'highs-ipm', 'highs',
            'interior-point', 'glpk' and 'cvxopt'.
            The defualt is 'ecos'.

        Returns
        -------
        pandas.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        return super().set_model(coef=coef, rtype=rtype, 
                                 mu=mu, mu0=mu0, aversion=aversion, ww0=ww0,
                                 hlength=hlength, method=method)


    def _wwgen(self):
        return MADAnalyzer(coef=self.coef, rtype=self.rtype,
                           method=self.method)
