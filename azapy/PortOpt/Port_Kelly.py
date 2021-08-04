# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 11:09:21 2021

@author: mircea
"""

from .Port_InvVol import Port_InvVol
from .KellyEngine import KellyEngine

class Port_Kelly(Port_InvVol):
    """
    Portfolio with Kelly optimal weights, periodically rebalanced.
    Functions: \n
        set_model \n
        get_port \n
        get_nshares \n
        get_weights \n
        get_account \n
        get_mktdata \n
        port_view \n
        port_view_all \n
        port_drawdown \n
        port_perf \n
        port_annual_returns \n
        port_monthly_returns
    """
    def set_model(self, rtype='Full', hlength=1.25):
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
            of months. The default is 1.25. 

        Returns
        -------
        pd.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        self._set_rtype(rtype)
        
        return super().set_model(hlength)
    
    def _set_rtype(self, rtype):
        valid_rtypes = ['Full', 'Order2']
        if rtype in valid_rtypes:
            self.rtype = rtype
        else:
            ValueError(f"Wrong rtype - must be one of {valid_rtypes}")
            
    def _ww_calc(self, data):
        return KellyEngine().getWeights(rrate=data, rtype=self.rtype)
