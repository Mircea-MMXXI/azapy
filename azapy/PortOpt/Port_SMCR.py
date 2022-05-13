from .Port_CVaR import Port_CVaR
from .SMCRAnalyzer import SMCRAnalyzer

class Port_SMCR(Port_CVaR):
    """
    Back testing the SMCR optimal portfolio periodically rebalanced.
    
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
    def set_model(self, mu, alpha=[0.9], coef=None, rtype='Sharpe',
                  hlength=3.25, method='ecos'):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        `mu` : float
            Reference rate. Its meaning depends on the value of rtype. For
            rtype equal to: \n
                'Sharpe' and 'Sharpe2': `mu` is the risk-free rate. \n
                'Ris' : `mu` is the targeted expected rate of returns. \n
                'MinRisk' and 'InvNrisk' : `mu` is ignored. \n
                'RiskAverse' : `mu` is the Lambda risk aversion coefficient.
        `alpha` : list, optional
            List of distinct alpha confidence levels. The default is [0.9].
        `coef` : list, optional
            List of positive mixture coefficients. Must have the same size with 
            `alpha`. A `None` value assumes an equal weighted risk mixture.
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
        `hlength` : float, optional
            The length in year of the historical calibration period relative
            to 'Dfix'. A fractional number will be rounded to an integer number
            of months. The default is 3.25 years.
        `method` : str, optional
            SOCP numerical method.
            It could be 'ecos' or 'cvxopt'.
            The default is 'ecos'.

        Returns
        -------
        pandas.DataFrame
            The portfolio time-series in the format 'date', 'pcolname'.
        """
        return super().set_model(mu, alpha, coef, rtype, hlength, method)
    
        
    def _wwgen(self):
        return SMCRAnalyzer(self.alpha, self.coef, rtype=self.rtype,
                            method=self.method)
