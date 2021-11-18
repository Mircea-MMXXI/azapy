import numpy as np
import pandas as pd

from .Port_ConstW import Port_ConstW

class Port_InvVol(Port_ConstW):
    """
    Back testing portfolio with weights proportional to the 
    inverse of component volatilities, periodically rebalanced.
    
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
    def __init__(self, mktdata, symb=None, sdate=None, edate=None, 
                 col_price='close', col_divd='divd', col_ref='adjusted',
                 col_calib='adjusted',
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
            should be present in mktdata. If set to None the symb will be 
            set to the full set of symbols present in mktdata. The default 
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
            information. The default is 'dvid'.
        col_ref : string, optional
            Column name in the mktdata DataFrame that will be used as a price 
            reference for portfolio components. The default is 'adjusted'.
        col_calib : string, optional
            Column name used for historical weights calibrations. 
            The default is 'adjusted'.
        pname : string, optional
            The name of the portfolio. The default is 'Port'.
        pcolname : string, optional
            Name of the portfolio price column. If it set to None than 
            pcolname=pname. The default is None.
        capital : float, optional
            Initial portfolio Capital in dollars. The default is 100000.
        schedule : pandas.DataFrame, optional
            Rebalancing schedule, with columns for 'Droll' rolling date and
            'Dfix' fixing date. If it is None than the schedule will be set 
            using the freq, nsoffset, fixoffset, hlength and calendar 
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
            business calendar. The default 
            is None.
    
        Returns
        -------
        The object.
        """
        super().__init__(mktdata=mktdata, symb=symb, 
                         sdate=sdate, edate=edate, 
                         col_price=col_price, col_divd=col_divd,
                         col_ref=col_ref, pname=pname,
                         pcolname=pcolname, capital=capital, 
                         schedule=schedule, freq=freq, noffset=noffset, 
                         fixoffset=fixoffset, calendar=calendar)
        self.col_calib = col_calib
    
    def set_model(self, hlength=3.25):
        """
        Set model parameters and evaluate the portfolio time-series.
        
        Parameters
        ----------
        hlength : float, optional
            The length in year of the historical calibration period relative 
            to 'Dfix'. A fractional number will be rounded to an integer number 
            of months. The default is 3.25 years. 

        Returns
        -------
        pd.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        self.hlength = hlength
        
        self._set_schedule()
        self._set_weights()
        self._port_calc()
        return self.port
    
    def _set_weights(self):
        mktdata = self.mktdata.pivot(columns='symbol', values=self.col_calib)
        periods = 63 if self.freq == 'Q' else 21
        
        # local function
        def _fww(rr):
            if rr.Dfix > self.edate:
                return pd.Series(np.nan, index=mktdata.columns)
            
            mm = mktdata[rr.Dhist:rr.Dfix].pct_change(periods=periods).dropna()
            return self._ww_calc(mm)
        
        w = self.schedule.apply(_fww, axis=1)
 
        self.ww = pd.concat([self.schedule, w], axis=1)
        
    def _ww_calc(self, data):
        vv = 1. / data.std()
        return vv / vv.sum()
