from .Port_InvVol import Port_InvVol
from .KellyEngine import KellyEngine

class Port_Kelly(Port_InvVol):
    """
    Backtesting Kelly optimal portfolio, periodically rebalanced.
    
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
    def set_model(self, rtype='ExpCone', hlength=1.25, method='ecos'):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        `rtype` : `str`, optional;
            Type of optimization. It could take the values:\n
                `'ExpCone'` - Exponential cone constraint programming solution
                for full Kelly problem. \n
                `'Full'` - Non-linear solver for full Kelly problem. \n
                `'Order2'` - Second order Tayler approximation of Kelly problem. \n
            The default is `'ExpCone'`.
        `hlength` : `float`, optional;
            The length in year of the historical calibration period relative 
            to `'Dfix'`. A fractional number will be rounded to an integer number 
            of months. The default is `1.25` years. 
        `method` : `str`, optional;
            The QP solver class. It is relevant only if `rtype='Order2'`.
            It takes 2 values: `'ecos'` or `'cvxopt'`.
            The default is `'ecos'`.

        Returns
        -------
        `pandas.DataFrame`;
            The portfolio time-series in the format "date", "pcolname".
        """
        self.rtype= rtype
        self.method = method
        
        return super().set_model(hlength)
 
            
    def _ww_calc(self, data):
        return KellyEngine().getWeights(rrate=data, 
                                        rtype=self.rtype, 
                                        method=self.method)
