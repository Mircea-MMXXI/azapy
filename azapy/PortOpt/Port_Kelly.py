from .Port_InvVol import Port_InvVol
from .KellyEngine import KellyEngine

class Port_Kelly(Port_InvVol):
    """
    Back testing portfolio with Kelly optimal weights, periodically rebalanced.
    
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
    def set_model(self, rtype='Full', hlength=1.25, method='ecos'):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        rtype : string, optional
            Type of optimization. It could take the values:\n
                "Full" - Non-linear (full) Kelly selection algorithm. \n
                "Order2" - Second order approximation of Kelly selection 
                algorithm. \n
            The default is 'Full'.
        hlength : float, optional
            The length in year of the historical calibration period relative 
            to 'Dfix'. A fractional number will be rounded to an integer number 
            of months. The default is 1.25 years. 
        method : string, optional
            The QP solver class. It is relevant only if rtype='Order2'.
            It takes 2 values: 'ecos' or None for default 'cvxopt' 
            algorithm.
            The default is 'ecos'.

        Returns
        -------
        pd.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        self._set_rtype(rtype)
        self.method = method
        
        return super().set_model(hlength)
    
    def _set_rtype(self, rtype):
        valid_rtypes = ['Full', 'Order2']
        if rtype in valid_rtypes:
            self.rtype = rtype
        else:
            ValueError(f"Wrong rtype - must be one of {valid_rtypes}")
            
    def _ww_calc(self, data):
        return KellyEngine().getWeights(rrate=data, rtype=self.rtype, 
                                        method=self.method)
