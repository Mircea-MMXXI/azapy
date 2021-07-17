# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 13:03:38 2021

@author: mircea
"""
import pandas as pd
import numpy as np

from .CVaRAnalyzer import CVaRAnalyzer
from .Port_InvVol import Port_InvVol


class Port_CVaR(Port_InvVol):
    """
    Portfolio with CVaR optimal weights, periodically rebalanced.
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
    def set_model(self, mu, alpha=[0.975], coef=None, rtype='Sharpe', 
                  hlength=3.25):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        mu : float
            Reference rate. Its meaning depends of the value of rtype. For
            rtype equal to: \n
                "Sharpe" : mu is the risk-free rate. \n
                "Risk" : mu is the targeted expected rate of returns. \n
                "MinRisk" and "InvNrisk" : mu is ignored. \n
                "RiskAverse" : mu is the Lambda risk aversion coefficient.
        alpha : list, optional
            The value of alpha CVaR confidence levels. The default is [0.975].
        coef : list, optional
            The coefficients values. The default is [1.].
        rtype : string, optional
            Type of optimization. It could take the values:\n
                "Sharpe" - Sharpe optimal portfolio. \n
                "Risk" - risk optimal portfolio. \n
                "MinRisk" - Minimum CVaR optimal portfolio. \n
                "InvNrisk" - optimal portfolio with same risk as the equally 
                weighted portfolio. \n
                "RiskAverse" - optimal portfolio for fixed risk aversion . \n
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
        self._set_alpha(alpha, coef)
        self._set_rtype(rtype)
        self.mu = mu
        self.hlength = hlength
        
        self._set_schedule()
        self._set_weights()
        self._port_calc()
        return self.port
    
    def _set_alpha(self, alpha, coef):
        # alpha
        self.alpha = np.array(alpha)
        assert np.all((0. < self.alpha) & (self.alpha < 1.)), \
            "alpha must be in (0, 1)"
        
        # coef
        if coef is None:
            self.coef = np.ones(len(self.alpha))
        else:
            assert len(coef) == len(self.alpha), \
                "coef must have same length as alpha"
            self.coef = np.array(coef)
            
        assert np.all(0. <= self.coef),  "coef must be >= 0"
        
        scoef = self.coef.sum()
        assert scoef > 0., "at leas one coef must be > 0"
        
        self.coef = self.coef / scoef
        
    def _set_rtype(self, rtype):
        rtype_values = ['Sharpe', 'Risk', 'MinRisk', 'InvNrisk', 'RiskAverse']
        assert rtype in rtype_values, \
            f"rtype must be one of {rtype_values}"
            
        self.rtype = rtype
        
    def _wwgen(self):
        return CVaRAnalyzer(self.alpha, self.coef, rtype=self.rtype)

    def _ww_calc(self, data):
        return self._wwgen().getWeights(mu=self.mu, rrate=data)