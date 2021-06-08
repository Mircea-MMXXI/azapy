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
                 pname='Omega_Port', pcolname=None, capital=100000, 
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
        
    def _set_wwgen(self):
        self.wwgen = OmegaAnalyzer(self.alpha[0], rtype=self.rtype)
        
    def get_port(self, mu=0., mu0=0., rtype='Sharpe'):
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
   
        self._make_schedule()
        self._make_ww()
        self._port_calc()
        return self.port
