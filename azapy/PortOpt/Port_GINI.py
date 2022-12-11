from .Port_CVaR import Port_CVaR
from .GINIAnalyzer import GINIAnalyzer

class Port_GINI(Port_CVaR):
    """
    Backtesting GINI optimal portfolio strategies, periodically rebalanced.

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
    def set_model(self, rtype='Sharpe', mu=None, mu0=0, aversion=None,
                  ww0=None, hlength=1.25, method='ecos', verbose=False):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        `rtype` : `str`, optional;
            Optimization type. Possible values: \n
                `'Risk'` : optimal-risk portfolio for targeted expected rate of 
                return.\n
                `'Sharpe'` : Sharpe-optimal portfolio - maximization solution.\n
                `'Sharpe2'` : Sharpe-optimal portfolio - minimization solution.\n
                `'MinRisk'` : minimum risk portfolio.\n
                `'RiskAverse'` : optimal-risk portfolio for a fixed 
                risk-aversion factor.\n
                `'InvNrisk'` : optimal-risk portfolio with the same risk value 
                as a benchmark portfolio (e.g., same as equal weighted 
                portfolio).\n
                `'Diverse'` : optimal-diversified portfolio for targeted
                expected rate of return (maximum of inverse 1-D).\n
                `'Diverse2'` : optimal-diversified portfolio for targeted
                expected rate of return (minmum of 1-D).\n
                `'MaxDiverse'` : maximum diversified portfolio.\n
                `'InvNdiverse'` : optimal-diversified portfolio with the same
                diversification factor as a benchmark portfolio 
                (e.g., same as equal weighted portfolio).\n
                `'InvNdrr'` : optima- diversified portfolio with the same 
                expected rate of return as a benchmark portfolio
                (e.g., same as equal weighted portfolio).\n
            The defauls is `'Sharpe'`.
        `mu` : `float`, optional;
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype='Risk'`
            The default is `None`.
        `mu0` : `float`, optional;
            Risk-free rate accessible to the investor.
            Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        `aversion` : `float`, optional;
            The value of the risk-aversion factor.
            Must be positive. Relevant only if `rtype='RiskAvers'`.
            The default is `None`.
        `ww0` : `list` (also `numpy.array` or `pandas.Series`), optional;
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
        `hlength` : `float`, optional;
            The length in year of the historical calibration period relative
            to `'Dfix'`. A fractional number will be rounded to an integer 
            number of months. The default is `3.25` years.
        `method` : `str`, optional;
            Linear programming numerical method.
            Could be: `'ecos'`, `'highs-ds'`, `'highs-ipm'`, `'highs'`,
            `'interior-point'`, `'glpk'` and `'cvxopt'`.
            The default is `'ecos'`.
        `verbose` : Boolean, optiona;
            If it set to `True` then it will print messages when the optimal
            portfolio degenerates to a single asset portfolio as a limited 
            case. 
            The default is `False`.

         Returns
        -------
        `pandas.DataFrame`;
            The portfolio time-series in the format 'date', 'pcolname'.
        """
        return super().set_model(rtype=rtype, mu=mu, mu0=mu0, 
                                 aversion=aversion, ww0=ww0, hlength=hlength, 
                                 method=method, verbose=verbose)


    def _wwgen(self):
        return GINIAnalyzer(rtype=self.rtype, method=self.method,
                            name=self.pname)
