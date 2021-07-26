# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 20:59:46 2021

@author: mircea
"""
from .Port_CVaR import Port_CVaR
from .OmegaAnalyzer import OmegaAnalyzer

class Port_Omega(Port_CVaR):
    """
    Portfolio with Omega optimal weights, periodically rebalanced.
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
    def set_model(self, mu=0., mu0=0., rtype='Sharpe', hlength=3.25,
                  method='ecos'):
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
        mu0 : float, optional
            Omega threshold rate (e.g. risk-free rate). The default is 0.
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
            Linear programming numerical method. 
            Could be one of 'ecos', 'highs-ds', 'highs-ipm', 'highs', 
            'interior-point', 'glpk' and 'cvxopt'.
            The defualt is 'ecos'.

        Returns
        -------
        pd.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        self._set_rtype(rtype)
        self.alpha = [mu0]
        self.coef = [1.]
        if self.rtype == 'Sharpe':
            self.mu = mu0
        else:
            self.mu = mu
        self.hlength = hlength
        self._set_method(method)
   
        self._set_schedule()
        self._set_weights()
        self._port_calc()
        return self.port

    def _wwgen(self):
        return OmegaAnalyzer(self.alpha[0], rtype=self.rtype, 
                             method=self.method)