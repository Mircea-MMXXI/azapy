# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 12:23:59 2021

@author: mircea
"""
import pandas as pd
import numpy as np

from .Port_Rebalanced import Port_Rebalanced
from azapy.MkT.readMkTData import NYSEgen
from azapy.util.schedule import schedule_roll

class Port_ConstW(Port_Rebalanced):
    """
    Portfolio with constant weights periodically rebalanced.
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
    def __init__(self, mktdata, symb=None, sdate=None, edate=None, 
                 col_price='close', col_divd='divd', col_ref='adjusted',
                 pname='Port', pcolname=None, capital=100000, 
                 schedule=None,
                 freq='Q', noffset=-3, fixoffset=-1, calendar=None):
        """
        Constructor

        Parameters
        ----------
        mktdata : pd.DataFrame
            MkT data in the format "symbol", "date", "open", "high", "low",
            "close", "volume", "adjusted", "divd", "split" (e.g. as returned
            by azapy.readMkT).
        symb : list, optional
            List of symbols for the basket components. All symbols MkT data
            should be included in mktdata. If set to None the symb will be 
            set to the full set of symbols included in mktdata. The default 
            is None.
        sdate : datetime, optional
            Start date for historical data. If set to None the sdate will 
            be set to the earliest date in mktdata. The default is None.
        edate : datetime, optional
            End date for historical dates and so the simulation. Must be 
            greater than  sdate. If it is None then edate will be set
            to the latest date in mktdata. The default is None.
        col_price : string, optional
            Column name in the mktdata DataFrame that will be considered 
            for portfolio aggregation.The default is 'close'.
        col_divd :  string, optional
            Column name in the mktdata DataFrame that holds the dividend 
            information. The default is 'dvid'
        col_ref : string, optional
            Column name in the mktdata DataFrame that will be used as a price 
            reference for portfolio components. The default is 'adjusted'.
        pname : string, optional
            The name of the portfolio. The default is 'Simple'.
        pcolname : string, optional
            Name of the portfolio price column. If it set to None that 
            pcolname=pname. The default is None.
        capital : float, optional
            Initial portfolio Capital in dollars. The default is 100000.
        schedule : pandas.DataFrame, optional
            Rebalancing schedule, with columns for 'Droll' rolling date and
            'Dfix' fixing date. If it is None than the schedule will be set 
            using the freq, nsoffset, fixoffset and calendar 
            information. The default is None.
        freq : string, optional
            rebalancing frequency. It can be 'Q' for quarterly or 'M' for 
            monthly rebalancing, respectively. It is relevant only is schedule 
            is None. The default is 'Q'.
        noffset : int, optional
            Number of business days offset for rebalancing date 'Droll' 
            relative to the end of the period (quart or month). A positive
            value add business days beyond the calendar end of the period while
            a negative value subtract business days. It is relevant only is 
            schedule is None. The default is -3.
        fixoffset : int, optional
            Number of business day offset of fixing date 'Dfix' relative to 
            the rebalancing date 'Droll'. It cane be 0 or negative. It is 
            relevant only is schedule is None. The default is -1.
        calendar : numpy.busdaycalendar, optional
            Business calendar. If it is None then it will be set to NYSE 
            business calendar via azapy.NYSEgen() function. The default 
            vale is None.

        Returns
        -------
        The object.
        """
        super().__init__(mktdata=mktdata, symb=symb, 
                         sdate=sdate, edate=edate, 
                         col_price=col_price, col_divd=col_divd,
                         col_ref=col_ref,
                         pname=pname,
                         pcolname=pcolname, capital=capital)
        self.schedule = schedule
        self.freq = freq
        self.noffset = noffset
        self.fixoffset = fixoffset
        self.calendar =  calendar 
        if self.calendar is None: self._default_calendar()

        
    def set_model(self, ww=None):
        """
         Set model parameters and evaluate the portfolio time-series.

        Parameters
        ----------
        ww : list (numpy.array to pandas.Series), optional
            List of weights. If it is panda.Series the index should match 
            the basket symb. Otherwise the weights are considered in the symb 
            order. If it is set to None than ww will be set to equal weights.
            The default is None.

        Returns
        -------
        pd.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        self.hlength = 0.
        
        self._set_schedule()
        self._set_weights(ww)
        self._port_calc()
        return self.port
    
    def _default_calendar(self):
        self.calendar = NYSEgen()
        
    def _set_schedule(self):
        if self.schedule is None:
            self.schedule = schedule_roll(self.sdate, self.edate, self.freq,
                                      self.noffset, self.fixoffset, 
                                      self.calendar, self.hlength)
            
    def _set_weights(self, ww):
        if ww is None:
            _ww = pd.Series(1., index=self.symb)
        elif isinstance(ww, pd.core.series.Series):
            _ww = ww
        else:
            _ww = pd.Series(ww, index=self.symb)
            
        assert _ww.size == self.symb.size, \
            f"ww size must have = number of symbols {self.symb.size}"
        assert np.all(_ww >= 0.), \
            "ww elements must be >= 0"
        wws = _ww.sum()
        assert wws > 0, "at least one ww element must be > 0"
        
        self.ww = self.schedule
        for sy in self.symb:
            self.ww[sy] = _ww[sy] / wws
