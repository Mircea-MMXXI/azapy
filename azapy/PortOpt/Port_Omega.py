from .Port_CVaR import Port_CVaR
from .OmegaAnalyzer import OmegaAnalyzer

class Port_Omega(Port_CVaR):
    """
    Back testing the Omega optimal portfolio weights, periodically rebalanced.

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
    def set_model(self, mu, alpha0=0., rtype='Sharpe', hlength=3.25,
                  method='ecos'):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        mu : float
            Reference rate. Its meaning depends on the value of `rtype`. For
            `rtype` equal to: \n
                "Sharpe" : `mu` is the risk-free rate \n
                "Risk" : `mu` is the targeted expected rate of returns \n
                "MinRisk" and "InvNrisk" : `mu` is ignored \n
                "RiskAverse" : `mu` is the Lambda risk aversion coefficient.
        alpha0 : float, optional
            Omega threshold rate (e.g. risk-free rate). The default is 0.
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
            The defualt is 'ecos'.

        Returns
        -------
        pandas.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        self._set_rtype(rtype)
        self.alpha = [alpha0]
        self.coef = [1.]
        self.mu = mu
        self.hlength = hlength
        self._set_method(method)

        self._set_schedule()
        self._set_weights()
        self._port_calc()
        return self.port

    def _wwgen(self):
        return OmegaAnalyzer(self.alpha[0], rtype=self.rtype,
                             method=self.method)
