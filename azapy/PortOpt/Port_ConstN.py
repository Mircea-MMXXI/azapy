# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 12:23:59 2021

@author: mirce
"""
import pandas as pd
import numpy as np
import pandas.tseries.offsets as pt

from .Port_Weighted import Port_Weighted
from azapy.MkT.readMkTData import NYSEgen

class Port_ConstN(Port_Weighted):
    """
    Portfolio with constant weights periodicaly rebalanced.
    Iherited from azapy.Port_Weighted \n
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
                 pname='ConstN', pcolname=None, capital=100000, 
                 freq='Q', noffset=-3, calendar=None):
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
            larget than  sdate. If it iset to None then edate will be sat
            to the latest date in rprice. The default is None.
        col : string, optional
            Name of column in the rprice DataFrame that will be considered 
            for portfolio agregation.The default is 'close'.
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
        noffset : intE, optional
            Number of offset business day form the calander end of invetment 
            period (rebalancing period). A positive value will add business 
            days beyond the calendar end of the period while a negative value
            will subtract business days. The default is -3.
        calendar : numpy.busdaycalendar, optional
            Business calendar compatible with the MkT data from rprice. If it
            None then it will be set to NYSE bunsiness calendar.
            The default is None.

        Returns
        -------
        The object.
        """
        super().__init__(rprice=rprice, symb=symb, 
                         sdate=sdate, edate=edate, 
                         col=col, pname=pname,
                         pcolname=pcolname, capital=capital)
        self.freq = freq
        self.noffset = noffset
        self.calendar =  calendar 
        if self.calendar is None: self._default_calendar()

        
    def get_port(self, ww=None):
        """
        Evaluates the portfolio timeseries.

        Parameters
        ----------
        ww : list (numpy.array ot pandas.Series), optional
            List of weights. If it is panda.Series the index should match 
            the basket symb. Othrwise the weights are considered in the symb 
            order. If it is set to None than ww will be set to equal weights.
            The default is None.

        Returns
        -------
        pd.DataFrame
            The portfolio timeseries in the format "date", "pcolname".
        """
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
        
        self.ww = self._make_simple_schedule()

        for sy in self.symb:
            self.ww[sy] = _ww[sy] / wws
 
        self._port_calc()
        return self.port
    
    def _make_simple_schedule(self):
        if self.freq == 'Q': edate = self.edate + pt.QuarterEnd(1)
        elif self.freq == 'M': edate = self.edate + pt.MonthEnd(1)
        else: raise ValueError("Wrong freq, Must be 'Q' or 'M'")
        
        tedx = pd.date_range(start=self.sdate, end=edate, freq=self.freq)\
            .to_numpy(dtype='<M8[D]')
        troll = np.busday_offset(tedx, self.noffset, roll='backward', 
                                 busdaycal=self.calendar)
        tfix = np.busday_offset(troll, -1, roll='backward', 
                                busdaycal=self.calendar)
        return pd.DataFrame({'Droll': troll, 'Dfix': tfix})
    
    def _default_calendar(self):
        self.calendar = NYSEgen()
