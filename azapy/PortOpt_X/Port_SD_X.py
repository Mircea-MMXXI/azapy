from .Port_MV_X import Port_MV_X
from azapy.Analyzers.SDAnalyzer import SDAnalyzer

class Port_SD_X(Port_MV_X):
    """
    Backtesting the SD optimal portfolio periodically rebalanced.
    
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
        return SDAnalyzer(name=self.pname,
                          rtype=self.rtype, mu=self.mu, mu0=self.mu0,
                          aversion=self.aversion, ww0=self.ww0,
                          method=self.method)
    