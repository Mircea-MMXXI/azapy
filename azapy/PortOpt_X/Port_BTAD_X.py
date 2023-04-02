from .Port_CVaR_X import Port_CVaR_X
from azapy.Analyzers.BTADAnalyzer import BTADAnalyzer

class Port_BTAD_X(Port_CVaR_X):
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
                  detrended=True, hlength=3.25, method='ecos', 
                  verbose=False):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        `alpha` : `list`, optional;
            List of thresholds. The default is `[0.]`.
        `coef` : `list`, optional;
            List of positive mixture coefficients. Note that `len(coef)`
            must be equal to `len(alpha)`. A `None` value assumes an
            equal weighted risk mixture.
            The vector of coefficients will be normalized to unit.
            The default is `None`.
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
        `detrendent` : Boolean, optional;
            If it set to `True` then the rates of return are detrended 
            (mean=0). The default value is `True`. 
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
        self.detrended = detrended
        return super().set_model(alpha=alpha, coef=coef, rtype=rtype, 
                                 mu=mu, mu0=mu0, aversion=aversion, ww0=ww0,
                                 hlength=hlength, method=method,
                                 verbose=verbose)


    def _wwgen(self):
        # return BTADAnalyzer(self.alpha, self.coef, rtype=self.rtype,
        #                     detrended=self.detrended, method=self.method,
        #                     name=self.pname)
        return BTADAnalyzer(self.alpha, self.coef, name=self.pname,
                            rtype=self.rtype, mu=self.mu, mu0=self.mu0,
                            aversion=self.aversion, ww0=self.ww0,
                            detrended=self.detrended, method=self.method)