import numpy as np
import pandas as pd
from azapy.Generators.Port_Generator import Port_Generator


class Port_Rebalanced(Port_Generator):
    """
    Backtesting a portfolio periodically rebalanced 
    (with an external schedule of weights).
    
    **Attributes**
        * `pname` : `str` - portfolio name
        * `ww` : `pandasDataFrame` - portfolio weights at each rebalancing date
        * `port` : `pandas.Series` - portfolio historical time-series
        * `schedule` : `pandas.DataFrame` - rebalancing schedule
       
    The most important method is `set_model`. It must be called before any
    other method.
    """                  
    def __init__(self, mktdata, symb=None, sdate=None, edate=None, 
                 col_price='close', col_divd='divd', col_ref='adjusted',
                 pname='Port', pcolname=None, capital=100000, schedule=None,
                 multitreading=True):
        """
        Constructor
    
        Parameters
        ----------
        mktdata : `pandas.DataFrame`
            MkT data in the format "symbol", "date", "open", "high", "low",
            "close", "volume", "adjusted", "divd", "split" (e.g., as returned
            by `azapy.readMkT` function).
        symb : `list`, optional
            List of symbols for the basket components. All symbols MkT data
            should be included in mktdata. If set to `None` the `symb` will be
            set to include all the symbols from `mktdata`. The default
            is `None`.
        sdate : date like, optional
            Start date for historical data. If set to `None` the `sdate` will
            be set to the earliest date in mktdata. The default is `None`.
        edate : date like, optional
            End date for historical dates and so the simulation. Must be
            greater than  `sdate`. If it is `None` then `edat`e will be set
            to the latest date in mktdata. The default is `None`.
        col_price : `str`, optional
            Column name in the mktdata DataFrame that will be considered
            for portfolio aggregation. The default is `'close'`.
        col_divd :  `str`, optional
            Column name in the mktdata DataFrame that holds the dividend
            information. The default is `'dvid'`.
        col_ref : `str`, optional
            Column name in the mktdata DataFrame that will be used as a price
            reference for portfolio components. The default is `'adjusted'`.
        pname : `str`, optional
            The name of the portfolio. The default is `'Port'`.
        pcolname : `str`, optional
            Name of the portfolio price column. If it set to `None` then
            `pcolname=pname`. The default is `None`.
        capital : `float`, optional
            Initial portfolio Capital in dollars. The default is `100000`.
        schedule : `pandas.DataFrame`, optional
            Rebalancing schedule, with columns for `'Droll'` rolling date and
            `'Dfix'` fixing date. If it is `None` than the schedule will be set
            using the `freq`, `noffset`, `fixoffset` and `calendar`
            information. The default is `None`.
        multitreading : Boolean, optional
            If it is `True` then the weights at the rebalancing dates will 
            be computed concurrent. The default is `True`.
    
        Returns
        -------
        The object.
        """
        super().__init__(mktdata=mktdata, symb=symb, sdate=sdate, edate=edate, 
                         col_price=col_price, col_divd=col_divd, 
                         col_ref=col_ref, pname=pname, pcolname=pcolname, 
                         capital=capital, schedule=schedule)
        self.nshares = None
        self.cash_invst = None
        self.cash_roll = None
        self.cash_divd = None
        self.schedule = schedule
        self.verbose = False
        
        
    def set_model(self, schedule=None, verbose=False):
        """
        Sets model parameters and evaluates the portfolio time-series.
        
        Parameters
        ----------
        schedule : `pandas.DataFrame`, optional
            Rebalancing schedule, with columns for `'Droll'` rolling date and
            `'Dfix'` fixing date. If it is `None` than the schedule will be set
            using the `freq`, `noffset`, `fixoffset` and `calendar`
            information. It is set to `None` it will overwrite the value 
            set by the constructor. The default is `None`.
        
        verbose : Boolean, optional:
            Sets teh verbose mode.

        Returns
        -------
        `pandas.DataFrame` : The portfolio time-series in the format 'date', 
        'pcolname'.
        """
        if schedule is not None:
            self.schedule = schedule
            
        self.ww = self.schedule
        self.verbose = verbose
        self._port_calc()
        return self.port
    
    
    def _port_calc(self):
        def _rank(x, ww=self.ww):
            for k in range(len(ww)):
                if x < ww.Droll[k]: return k
            return len(ww)
        
        mktdata = self.mktdata.pivot(columns='symbol', values=self.col_price)
        div = self.mktdata.pivot(columns='symbol', values=self.col_divd)
        lw = np.zeros([div.shape[0]], dtype=int)
        for dx in self.ww.Droll:
            lw[div.index >= dx] += 1
        mmix = pd.MultiIndex.from_arrays([lw, div.index], names=('lw', 'date'))
        div.index = mmix
        div = div.groupby(level='lw').sum()
        symb = self.symb

        mktdata.index = mmix
        mktgr = mktdata.groupby(level='lw')
        mktdata = mktdata.droplevel(0)

        self.port = []
        self.nshares = []    
        self.cash_invst = []
        self.cash_roll = []
        self.cash_divd = [0.]   
        cap = self.capital
        for k, v in mktgr:
            if k == 0:  continue
            v = v.droplevel(0)
            
            nsh = (self.ww[symb].iloc[k - 1] * cap
                    / mktdata.loc[self.ww.Dfix[k - 1]]).round(0)
            self.nshares.append(nsh)
            self.port.append(v @ nsh)
            
            invst = nsh @ mktdata.loc[self.ww.Droll[k - 1]]
            dcap = cap - invst
            divd = div.loc[k] @ nsh
            cap = self.port[-1][-1] + divd  + dcap
            
            self.cash_invst.append(invst)
            self.cash_roll.append(dcap)
            self.cash_divd.append(divd)

        self.port = pd.concat(self.port) \
            .pipe(pd.DataFrame, columns=[self.pcolname])

        self.nshares = pd.DataFrame(self.nshares,
                                    index=self.ww.Droll[:len(self.nshares)]) 
    