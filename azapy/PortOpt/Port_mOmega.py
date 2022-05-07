import numpy as np
from .Port_CVaR import Port_CVaR
from .mOmegaAnalyzer import mOmegaAnalyzer

class Port_mOmega(Port_CVaR):
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
    def set_model(self, mu, alpha=[0.], coef=[1.], rtype='Sharpe', 
                  hlength=3.25, method='ecos'):
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
        alpha : list, optional
            Omega threshold levels. The default is [0.].
        coef : list, optional
            List of mixture coefficients values. Note that `len(coef)`
            must be equal to `len(alpha)`. A value of `None` assumes
            equal weights,
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
            The defualt is 'ecos'.

        Returns
        -------
        pandas.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        return super().set_model(mu, alpha, coef, rtype, hlength, method)

    
    def _set_alpha(self, alpha, coef):
        # alpha
        if len(alpha) <= 0:
            raise ValueError("alpha is an empty list")
            
        self.alpha = np.array(alpha)

        # coef
        if coef is None:
            self.coef = np.ones(len(self.alpha))
        else:
            if len(coef) != len(self.alpha):
                raise ValueError("coef must have same length as alpha")
            self.coef = np.array(coef)

        if np.any(self.coef <= 0.):
            raise ValueError("coef must be > 0")

        scoef = self.coef.sum()
        self.coef = self.coef / scoef


    def _wwgen(self):
        return mOmegaAnalyzer(self.alpha, self.coef, rtype=self.rtype,
                             method=self.method)
