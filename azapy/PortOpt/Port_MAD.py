# -*- coding: utf-8 -*-
"""
Created on Mon May 17 11:12:37 2021

@author: mircea
"""
import numpy as np

from .Port_CVaR import Port_CVaR
from .MADAnalyzer import MADAnalyzer

class Port_MAD(Port_CVaR):
    """
    Portfolio with MAD optimal weights, periodically rebalanced.
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
    def set_model(self, mu, coef=[1.], rtype='Sharpe', hlength=3.25):
        """
        Set model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        mu : float
            Reference rate. Its meaning depends of the value of rtype. For
            rtype equal to: \n
                "Sharpe" : mu is the risk-free rate \n
                "Risk" : mu is the targeted expected rate of returns \n
                "MinRisk" and "InvNrisk" : mu is ignored
        coef : list, optional
            The coefficients values. The default is [1.].
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

        Returns
        -------
        pd.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        return super().set_model(mu=mu, coef=coef, rtype=rtype, 
                                 hlength=hlength)
    
    def _set_alpha(self, alpha, coef):
        # ignore alpha
        coef = np.trim_zeros(np.array(coef), trim='b')
        assert np.all(coef >= 0.), "all coef must be non-negative"
        ssc = np.sum(coef)
        assert ssc > 0., "at least one coef must be > 0"
        self.coef = coef / ssc
        
    def _set_wwgen(self):
        self.wwgen = MADAnalyzer(coef=self.coef, rtype=self.rtype)
