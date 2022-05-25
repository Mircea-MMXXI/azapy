from .Port_CVaR import Port_CVaR
from .BTADAnalyzer import BTADAnalyzer

class Port_BTAD(Port_CVaR):
    """
    Backtesting the BTAD optimal portfolio strategies, periodically rebalanced.

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
    def set_model(self, mu, alpha=[0.], coef=None, rtype='Sharpe', 
                  detrended=False, hlength=3.25, method='ecos'):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        `mu` : float
            Reference rate. Its meaning depends on the value of `rtype`. For
            `rtype` equal to: \n
                "Sharpe" and "Sharpe2": `mu` is the risk-free rate \n
                "Risk" : `mu` is the targeted expected rate of returns \n
                "MinRisk" and "InvNrisk" : `mu` is ignored \n
                "RiskAverse" : `mu` is the Lambda risk aversion coefficient.
        `alpha` : list, optional
            List of Omega thresholds. The default is [0.].
        `coef` : list, optional
            List of positive mixture 
            coefficients. Must have the same size as `alpha`. 
            A `None` value assumes an equal weighted risk mixture.
            The vector of coefficients will be normalized to unit.
            The default is `None`.
        `rtype` : str, optional
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
        `detrended` : Boolean, optional
            In the Delta-risk expression use: \n
                `True` : detrended rate of return, i.e. r - E(r), \n
                `False` : standard rate of return. 
            The default is `False`.
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
        self.detrended = detrended
        return super().set_model(mu, alpha, coef, rtype, hlength, method)


    def _wwgen(self):
        return BTADAnalyzer(self.alpha, self.coef, rtype=self.rtype,
                            detrended=self.detrended, method=self.method)
