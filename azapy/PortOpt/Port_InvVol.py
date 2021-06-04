# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 16:11:18 2021

@author: mircea
"""
import pandas as pd
import numpy as np

from .Port_ConstN import Port_ConstN

class Port_InvVol(Port_ConstN):
    """
    Portfolio with weights proportional to inverse of volatility, 
    periodicaly rebalanced.
    Iherited from azapy.Port_ConstN \n
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
                 pname='InvVol', pcolname=None, capital=100000, 
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
        noffset : int, optional
            Number of offset business day form the calander end of invetment 
            period (rebalancing period). A positive value will add business 
            days beyond the calendar end of the period while a negative value
            will subtract business days. The default is -3.
        hlenght : float, optional
            Defines the calibration period in years for basket component 
            volatilities. The calibration period is prior and ends on the 
            fixing date. It could be a fractionar number but the actual 
            calibration period will rounded to the nearest multiple of 
            rebalancing periods. The defualt is 1.
        calendar : numpy.busdaycalendar, optional
            Business calendar compatible with the MkT data from rprice. If it
            None then it will be set to NYSE bunsiness calendar.
            The default is None.

        Returns
        -------
        The object.
        """
        super().__init__(rprice, symb, sdate, edate, col, 
                        pname, pcolname, capital, freq, noffset, calendar)
        self.hlenght = hlenght
        
        self.schedule = None
        
    def get_port(self):
        """
        Evaluates the portfolio timeseries.

        Returns
        -------
        pd.DataFrame
            The portfolio timeseries in the format "date", "pcolname".
        """
        self._make_schedule()
        self._make_ww()
        self._port_calc()
        return self.port
    
    def _make_schedule(self):
        sch = self._make_simple_schedule()
        sch['Dhist'] = sch['Dfix'] - pd.offsets \
            .DateOffset(months=round(self.hlenght * 12, 0))
        sch['Dhist'] = np.busday_offset(sch['Dhist'].to_numpy(dtype='<M8[D]'), 
                                        0, roll='backward', 
                                        busdaycal=self.calendar)
        for k in range(len(sch)):
            if sch['Dhist'][k] >= self.sdate: 
                self.schedule = sch[k:].reset_index(drop=True)
                return 
        
        raise ValueError("Cannot make a schedule!!")
        
    def _make_ww(self):
        mktdata = self.mktdata.pivot(columns='symbol', values='adjusted')
            
        w = []
        for isch in self.schedule.index:
            md = mktdata[self.schedule.Dhist[isch]:self.schedule.Droll[isch]]
            vv = 1. / md.diff().std()
            w.append(vv / vv.sum())
        
        self.ww = pd.concat([self.schedule, pd.DataFrame(w)], axis=1)
        print(self.ww)
