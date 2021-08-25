from .Port_MAD import Port_MAD
from .LSSDAnalyzer import LSSDAnalyzer

class Port_LSSD(Port_MAD):
    """
    Backtesting the LSSD optimal portfolio periodically rebalanced.
    
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
    def _set_method(self, method):
        methods = ['ecos', 'cvxopt']
        if not method in methods:
            raise ValueError(f"mehtod must be one of {methods}")   
        self.method = method
        
        
    def _wwgen(self):
        return LSSDAnalyzer(coef=self.coef, rtype=self.rtype, 
                            method=self.method)
    