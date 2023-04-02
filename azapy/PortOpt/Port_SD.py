from .Port_MV import Port_MV
#from .SDAnalyzer import SDAnalyzer
from azapy.Analyzers.SDAnalyzer import SDAnalyzer

class Port_SD(Port_MV):
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
        return SDAnalyzer(freq=self.freq, 
                          hlength=self.hlength, calendar=self.calendar,
                          name=self.pname,
                          rtype=self.rtype, mu=self.mu, mu0=self.mu0,
                          aversion=self.aversion, ww0=self.ww0,
                          method=self.method)
    