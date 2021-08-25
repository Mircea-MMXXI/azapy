from .Port_CVaR import Port_CVaR
from .MVAnalyzer import MVAnalyzer

class Port_MV(Port_CVaR):
    """
    Backtesting the MV optimal portfolio periodically rebalanced.
    
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
    def set_model(self, mu, rtype='Sharpe', hlength=3.25, method='ecos'):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        mu : float
            Reference rate. Its meaning depends of the value of rtype. For
            rtype equal to: \n
                "Sharpe" : mu is the risk-free rate \n
                "Risk" : mu is the targeted expected rate of returns \n
                "MinRisk" and "InvNrisk" : mu is ignored
        rtype : string, optional
            Type of optimization. It could take the values:\n
                "Sharpe" - C-Sharpe optimal portfolio \n
                "Risk" - CVaR optimal portfolio \n
                "MinRisk" - Minimum CVaR optimal portfolio \n
                "InvNrisk" - optimal portfolio with same CVaR as the equally 
                weighted portfolio. \n
                The default is 'Sharpe'.
        hlength : float, optional
            The length in year of the historical calibration period relative 
            to 'Dfix'. A fractional number will be rounded to an integer number 
            of months. The default is 3.25. 
        method : string, optional
            Numerical method to solve the SOCP and QP. Can take one of the 
            values 'ecos' and 'cvxopt'. The default is 'ecos'.

        Returns
        -------
        pd.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        return super().set_model(mu=mu, rtype=rtype, hlength=hlength, 
                                 method=method)
    
    def _set_method(self, method):
        methods = ['ecos', 'cvxopt']
        if not method in methods:
            raise ValueError(f"mehtod must be one of {methods}")          
        self.method = method
        
    def _wwgen(self):
        return MVAnalyzer(rtype=self.rtype, method=self.method)
