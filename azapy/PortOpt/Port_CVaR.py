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
    Inherits from azapy.Port_InvVol \n
    Functions: \n
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
    def __init__(self, rprice, symb=None, sdate=None, edate=None, col='close', 
                 pname='CVaR_Port', pcolname=None, capital=100000, 
                 freq='Q', noffset=-3, hlenght=1, calendar=None):
        """
        Constructor

        Parameters
        ----------
        rprice : pd.DataFrame
            MkT data in the format "symbol", "date", "open", "high", "low",
            "close", "volume", "adjusted", "divd", "split" (e.g. as returned
            by azapy.readMkT).
        symb : list, optional
            List of symbols for the basket components. All symbols MkT data
            should be included in rprice. If set to None the symb will be 
            set to the full set of symbols included in rprice. The default 
            is None.
        sdate : datetime, optional
            Start date for historical data. If set to None the sdate will 
            be set to the earliest date in rprice. The default is None.
        edate : datetime, optional
            End date for historical dates and so the simulation. Must be 
            greater than  sdate. If it is None then edate will be set
            to the latest date in rprice. The default is None.
        col : string, optional
            Name of column in the rprice DataFrame that will be considered 
            for portfolio aggregation.The default is 'close'.
        pname : string, optional
            The name of the portfolio. The default is 'Simple'.
        pcolname : string, optional
            Name of the portfolio price column. If it set to None that 
            pcolname=pname. The default is None.
        capital : float, optional
            Initial portfolio Capital in dollars. The default is 100000.
        freq : string, optional
            Defines the rebalancing period. Can take the following values:
                "M" : monthly rebalancing \n
                "Q" : quarterly rebalancing \n
                The default is 'Q'. 
        noffset : int, optional
            Number of offset business day form the calendar end of investment 
            period (rebalancing period). A positive value will add business 
            days beyond the calendar end of the period while a negative value
            will subtract business days. The default is -3.
        hlenght : float, optional
            Defines the calibration period in years for basket component 
            volatilities. The calibration period is prior and ends on the 
            fixing date. It could be a fractional number but the actual 
            calibration period will rounded to the nearest multiple of 
            rebalancing periods. The default is 1.
        calendar : numpy.busdaycalendar, optional
            Business calendar compatible with the MkT data from rprice. If it
            None then it will be set to NYSE business calendar.
            The default is None.

        Returns
        -------
        The object.
        """
        super().__init__(rprice, symb, sdate, edate, col, 
                         pname, pcolname, capital, freq, noffset, hlenght,
                         calendar)
        self.mu = None
        self.alpha = None
        self.coef = None
        self.rtype = None
        self.wwgen = None
        
    def get_port(self, mu, alpha=[0.975], coef=[1.], rtype='Sharpe'):
        """
        Evaluates the portfolio timeseries.

        Parameters
        ----------
        mu : float
            Reference rate. Its meaning depends of the value of rtype. For
            rtype equal to: \n
                "Sharpe" : mu is the risk-free rate \n
                "Risk" : mu is the targeted expected rate of returns \n
                "MinRisk" and "InvNrisk" : mu is ignored
        alpha : list, optional
            The value of alpha CVaR confidence levels. The default is [0.975].
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

        Returns
        -------
        pd.DataFrame
            The portfolio time-series in the format "date", "pcolname".

        """
        self._set_alpha(alpha, coef)
        self._set_rtype(rtype)
        self.mu = mu
        
        self._make_schedule()
        self._make_ww()
        self._port_calc()
        return self.port
    
    def _set_alpha(self, alpha, coef):
        self.alpha = np.array(alpha)
        assert np.all((0. < self.alpha) & (self.alpha < 1.)), \
            "alpha must be in (0, 1)"
        self.coef = np.array(coef)
        assert np.all(0. <= self.coef), \
            "coef must be >= 0"
        scoef = self.coef.sum()
        assert scoef > 0., "at leas one coef must be > 0"
        self.coef = self.coef / scoef
        
    def _set_rtype(self, rtype):
        rtype_values = ['Sharpe', 'Risk', 'MinRisk', 'InvNrisk']
        assert rtype in rtype_values, \
            f"rtype must be one of {rtype_values}"
        self.rtype = rtype
        
    def _set_wwgen(self):
        self.wwgen = CVaRAnalyzer(self.alpha, self.coef, rtype=self.rtype)
        
    def _make_ww(self):
        mktdata = self.mktdata.pivot(columns='symbol', values='adjusted')
        periods = 62 if self.freq == 'Q' else 21
        self._set_wwgen()
        
        w = []
        for isch in self.schedule.index:
            if self.schedule.Dfix[isch] > self.edate: 
                w.append(pd.Series(np.nan, index=self.mktdata.columns))
                continue
            md = mktdata[self.schedule.Dhist[isch]:self.schedule.Droll[isch]] \
                .pct_change(periods=periods).dropna()
            w.append(self.wwgen.getWeights(mu=self.mu, rrate=md))
        
        self.ww = pd.concat([self.schedule, 
                             pd.DataFrame(w, columns=mktdata.columns)], axis=1)
