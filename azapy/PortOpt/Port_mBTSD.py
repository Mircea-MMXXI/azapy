from .Port_mOmega import Port_mOmega
from .mBTSDAnalyzer import mBTSDAnalyzer

class Port_mBTSD(Port_mOmega):
    """
    Back testing the BTSD optimal portfolio weights, periodically rebalanced.

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
            BTSD threshold levels. The default is [0.].
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
            SOCP numerical method.
            Could be: 'ecos' or 'cvxopt'.
            The defualt is 'ecos'.

        Returns
        -------
        pandas.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        return super().set_model(mu, alpha, coef, rtype, hlength, method)
 
    
    def _wwgen(self):
        return mBTSDAnalyzer(self.alpha, self.coef, rtype=self.rtype,
                            method=self.method)