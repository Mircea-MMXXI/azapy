import numpy as np
import pandas as pd
import threading
import copy

from .Port_Simple import Port_Simple
from azapy.util.schedule import schedule_offset
from azapy.MkT.MkTcalendar import NYSEgen
from azapy.util.drawdown import max_drawdown


class Port_Generator(Port_Simple):
    """
    Backtesting a portfolio periodically rebalanced with model genrated 
    weights.
    
    Methods:
        * set_model
        * get_port
        * get_weights
        * get_nshares
        * get_account
        * get_mktdata
        * port_view
        * port_view_all
        * port_drawdown
        * port_perf
        * port_annual_returns
        * port_monthly_returns
        * port_period_returns
        * port_period_perf
    Attributs:
        * pname
        * ww
        * port
        * schedule
    """
    def __init__(self, mktdata, symb=None, sdate=None, edate=None, 
                 col_price='close', col_divd='divd', col_ref='adjusted',
                 pname='Port', pcolname=None, capital=100000, schedule=None,
                 freq='Q', noffset=-3, fixoffset=-1, histoffset=3.25, 
                 calendar=None, multitreading=True):
        """
        Constructor
    
        Parameters
        ----------
        `mktdata` : `pandas.DataFrame`;
            MkT data in the format "symbol", "date", "open", "high", "low",
            "close", "volume", "adjusted", "divd", "split" (e.g. as returned
            by `azapy.readMkT` function).
        `symb` : `list`, optional;
            List of symbols for the basket components. All symbols MkT data
            should be included in mktdata. If set to `None` the `symb` will be
            set to include all the symbols from `mktdata`. The default
            is `None`.
        `sdate` : date like, optional;
            Start date for historical data. If set to `None` the `sdate` will
            be set to the earliest date in mktdata. The default is `None`.
        `edate` : date like, optional;
            End date for historical dates and so the simulation. Must be
            greater than  `sdate`. If it is `None` then `edat`e will be set
            to the latest date in mktdata. The default is `None`.
        `col_price` : `str`, optional;
            Column name in the mktdata DataFrame that will be considered
            for portfolio aggregation. The default is `'close'`.
        `col_divd` :  `str`, optional;
            Column name in the mktdata DataFrame that holds the dividend
            information. The default is `'dvid'`.
        `col_ref` : `str`, optional;
            Column name in the mktdata DataFrame that will be used as a price
            reference for portfolio components. The default is `'adjusted'`.
        `pname` : `str`, optional;
            The name of the portfolio. The default is `'Port'`.
        `pcolname` : `str`, optional;
            Name of the portfolio price column. If it set to `None` then
            `pcolname=pname`. The default is `None`.
        `capital` : `float`, optional;
            Initial portfolio Capital in dollars. The default is `100000`.
        `schedule` : `pandas.DataFrame`, optional;
            Rebalancing schedule, with columns for `'Droll'` rolling date and
            `'Dfix'` fixing date. If it is `None` than the schedule will be set
            using the `freq`, `noffset`, `fixoffset` and `calendar`
            information. The default is `None`.
        `freq` : `str`, optional;
            rebalancing frequency. It can be `'Q'` for quarterly or `'M'` for
            monthly rebalancing, respectively. It is relevant only is schedule
            is `None`. The default is `'Q'`.
        `noffset` : `int`, optional;
            Number of business days offset for rebalancing date `'Droll'`
            relative to the end of the period (quart or month). A positive
            value add business days beyond the calendar end of the period while
            a negative value subtracts business days. It is relevant only is
            schedule is `None`. The default is `-3`.
        `fixoffset` : `int`, optional;
            Number of business day offset of fixing date `'Dfix'` relative to
            the rebalancing date `'Droll'`. It can be 0 or negative. It is
            relevant only is schedule is `None`. The default is `-1`.
        `calendar` : `numpy.busdaycalendar`, optional;
            Business calendar. If it is `None` then it will be set to NYSE
            business calendar. The default
            vale is `None`.
        `multitreading` : Boolean, optional;
            If it is `True` then the weights at the rebalancing dates will 
            be computed concurrent. The default is `True`.
    
        Returns
        -------
        The object.
        """
        super().__init__(mktdata=mktdata, symb=symb, 
                         sdate=sdate, edate=edate, 
                         col=col_ref, pname=pname, 
                         pcolname=pcolname, capital=capital)
        self.col_price = col_price
        self.col_divd = col_divd
        self.nshares = None
        self.cash_invst = None
        self.cash_roll = None
        self.cash_divd = None
        self.schedule = schedule
        self.freq = freq
        self.noffset = noffset
        self.fixoffset = fixoffset
        self.histoffset = histoffset
        self.calendar = calendar
        self.verbose = False
        if self.calendar is None: 
            self._default_calendar()
        self.multitreading = multitreading
        
        
    def set_model(self, wwModel, verbose=False):
        """
        Sets model parameters and evaluates the portfolio time-series.
        
        Parameters
        ----------
        `wwModel` : `ModelPipeline` object;
            Weights model
        `verbose` : Boolean, optional:
            Sets teh verbose mode.

        Returns
        -------
        `pandas.DataFrame`;
            The portfolio time-series in the format "date", "pcolname".
        """
        self.wwModel = wwModel
        self.verbose = verbose
        
        self._set_schedule()
        self._set_weights()
        self._port_calc()

        return self.port
    
    
    def _set_schedule(self):
        if self.schedule is None:
            self.schedule = schedule_offset(self.sdate, self.edate, self.freq,
                                            self.noffset, self.fixoffset,
                                            self.calendar, self.histoffset)
    
    
    def _set_weights(self):
        if self.multitreading:
            return self._set_weights_mt()
        return self._set_weights_sr()
    
    
    def _set_weights_sr(self):
        # local function
        def _fww(Dfix):
            if Dfix > self.edate:
                return pd.Series(np.nan, index=self.mktdata['symbol'].unique())
            
            mm = self.mktdata.loc[self.mktdata.index <= Dfix]
            return self.wwModel.getWeights(mm, pclose=False, 
                                           verbose=self.verbose)
        
        w = pd.DataFrame([_fww(Dfix) for Dfix in self.schedule['Dfix'].tolist()])
        self.ww = pd.concat([self.schedule, w], axis=1)


    def _set_weights_mt(self):
        # local function
        def _fww(Dfix, mod):
            if Dfix > self.edate:
                return pd.Series(np.nan, index=self.mktdata['symbol'].unique())
            
            mm = self.mktdata.loc[self.mktdata.index <= Dfix].copy()
            
            wei = mod.getWeights(mm, pclose=False, verbose=self.verbose)
            wei['Dfix'] = Dfix
            w.append(wei)
        
        w = []
        th = [threading.Thread(target=_fww, 
                               args=(k, copy.deepcopy(self.wwModel))) \
              for k in self.schedule['Dfix'].tolist()]
        for i in range(len(th)):
            th[i].start()
        for i in range(len(th)):
            th[i].join()
            
        self.ww = pd.merge(left=self.schedule, right=pd.DataFrame(w), on='Dfix')


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
        if '_CASH_' not in mktdata.columns:            
            mktdata['_CASH_'] = 1 
            div['_CASH_'] = 0 
            symb = np.append(symb, '_CASH_')
       
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
  
    
    def _default_calendar(self):
        self.calendar = NYSEgen()
        
    
    def get_account(self, fancy=False):
        """
        Returns additional bookkeeping information regarding rebalancing 
        (e.g. residual cash due rounding number of shares, previous period 
         dividend cash accumulation, etc.)

        Parameters
        ----------
        `fancy` : Boolean, optional;
            * `False`: the values are reported in unaltered algebraic format. 
            * `True` : the values are reported rounded.
            
        The default is `False`.

        Returns
        -------
        `pandas.DataFrame`;
            Reports, for each rolling period identified by `'Droll'`: 

            * number of shares hold for each symbol,
            * 'cash_invst' : cash invested at the beginning of period,
            * 'cash_roll' : cash rolled to the next period,
            * 'cash_divd' : cash dividend accumulated in the previous period.
                
        Note: The capital at the beginning of the period is 
        cash_invst + cash_roll. It is also equal to the previous period: 
        value of the shares on the fixing date + cash_roll + cash_divd.
        There are 2 sources for the cash_roll. The roundup to integer 
        number of shares and the shares close price differences between 
        the fixing (computation) and rolling (execution) dates. It could
        be positive or negative. The finance of the cash_roll during 
        each rolling period is assumed to be done separately by the investor.
        """
        acc_tab = self.nshares.copy()
        acc_tab['cash_invst'] = self.cash_invst
        acc_tab['cash_roll'] = self.cash_roll
        acc_tab['cash_divd'] = self.cash_divd[:-1]
        
        if not fancy: return acc_tab
        
        acc_tab = acc_tab.round(2)
        acc_tab[self.symb] = acc_tab[self.symb].astype('int')
        
        return acc_tab


    def port_period_returns(self, fancy=False):
        """
        Computes the rolling periods rate of returns. 

        Parameters
        ----------
        `fancy` : Boolean, optional;
           * `False`: returns in algebraic form.
           * `True`: returns percentage rounded to 2 decimals.
        The default is `False`.

        Returns
        -------
        `pandas.DataFrame`;
            Each rolling period is indicated by its start date, `'Droll'`. 
            Included are the fixing data, `'Dfix'`, and the 
            portfolio weights.
        """
        # local function
        def frr(x):
            p2 = x.p1.shift(-1)
            p2.iloc[-1] = self.port.iloc[-1,-1]
            return p2 / x.p1 - 1
        
        rww1 = self.ww.merge(self.port, left_on='Droll', right_index=True)
        rww1['Rx'] = rww1[self.port.columns[0]].pct_change()
        print(rww1)
        
        rww = self.ww.merge(self.port, left_on='Droll', right_index=True)\
                  .rename(columns={self.port.columns[0]: 'p1'})\
                  .assign(RR=frr)[['Droll', 'Dfix', 'RR'] + list(self.symb)]
 
        if not fancy:
            return rww
        
        rww['RR'] = rww['RR'].round(4) * 100
        rww[self.symb] = rww[self.symb].round(4).abs() * 100
        
        return rww
    
    
    def port_period_perf(self, fancy=False):
        """
        Returns portfolio performance for each rolling period i.e. 
        the rate of return, the rolling 
        min and max returns, and max drawdown during the period.

        Parameters
        ----------
        `fancy` : Boolean, optional;
           * `False`: returns in algebraic form.
           * `True`: returns percentage rounded to 2 decimals.
        The default is `False`.

        Returns
        -------
        `pandas.DataFrame`;
            Each rolling period is indicated by its start date, `'Droll'`.\n
            `'RR'` - period rate of return.\n
            `'RR_Min'` - minimum rolling rate of return in the period.\n
            `'RR_Max'` - maximum rolling rate of return in the period.\n
            `'DD_Max'` - maximum drawdown in the period.\n
            `'RR_Min_Date'` - date of `'RR_Min'`.\n
            `'RR_Max_Date'` - date of `'RR_Max'`.\n
            `'DD_Max_Date'` - date of `'DD_Max'`.\n
        """
        nn = self.schedule.shape[0]
        rr = np.zeros(nn)
        rmin = np.zeros(nn)
        rmin_date = []
        rmax = np.zeros(nn)
        rmax_date = []
        ddmax = np.zeros(nn)
        ddmax_date = []
        for i in range(nn):
            sd = self.schedule['Droll'][i]
            if i == (nn - 1):
                ed = self.port.index[-1]
            else:
                ed = self.schedule['Droll'][i + 1]

            pr = self.port.loc[sd:ed] / self.port.loc[sd] - 1
            rr[i] = pr.iloc[-1]
            mm = pr.idxmin(numeric_only=True)[0]
            rmin[i] = pr.loc[mm]
            rmin_date.append(mm)
            mm = pr.idxmax(numeric_only=True)[0]
            rmax[i] = pr.loc[mm]
            rmax_date.append(mm)
            dd = max_drawdown(self.port.loc[sd:ed], self.port.columns[0])
            ddmax[i] = dd[0]
            ddmax_date.append(dd[1])
            
        rout = pd.DataFrame({'Droll': self.schedule['Droll'],
                             'RR': rr, 
                             'RR_Min': rmin, 'RR_Max': rmax, 'DD_Max': ddmax,
                             'RR_Min_Date': rmin_date,
                             'RR_Max_Date': rmax_date,
                             'RR_DD_Date': ddmax_date
                             })
        
        if not fancy:
            return rout
        
        nll = ['RR', 'RR_Min', 'RR_Max', 'DD_Max']
        rout[nll] = rout[nll].round(4) * 100
        return rout
            