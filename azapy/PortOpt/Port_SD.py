from .Port_MV import Port_MV
from .SDAnalyzer import SDAnalyzer

class Port_SD(Port_MV):
    """
    Back testing the SD optimal portfolio periodically rebalanced.
    
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
    def _wwgen(self):
        return SDAnalyzer(rtype=self.rtype, method=self.method)
    