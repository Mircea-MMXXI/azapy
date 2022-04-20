import pandas as pd

from .Port_Simple import Port_Simple

class Port_Rebalanced(Port_Simple):
    """
    Back testing the portfolio with custom scheduled weights.
    
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
                 pname='Port', pcolname=None, capital=100000):
        """
        Constructor

        Parameters
        ----------
        mktdata : pandas.DataFrame
            MkT data in the format "symbol", "date", "open", "high", "low",
            "close", "volume", "adjusted", "divd", "split" (e.g. as returned
            by `azapy.readMkT` function).
        symb : list, optional
            List of symbols for the basket components. All symbols MkT data
            should be included in `mktdata`. If set to `None` the `symb` will be 
            set to the full set of symbols included in `mktdata`. The default 
            is `None`.
        sdate : date like, optional
            Start date for historical data. If set to `None` then the `sdate` will 
            be set to the earliest date in mktdata. The default is `None`.
        edate : date like, optional
            End date for historical dates and so the simulation. Must be 
            greater than  `sdate`. If it is `None` then `edate` will be set
            to the latest date in `mktdata`. The default is `None`.
        col_price : str, optional
            Column name in the mktdata DataFrame that will be considered 
            for portfolio aggregation. The default is 'close'.
        col_divd :  str, optional
            Column name in the mktdata DataFrame that holds the dividend 
            information. The default is 'dvid'.
        col_ref : str, optional
            Column name in the mktdata DataFrame that will be used as a price 
            reference for portfolio components. The default is 'adjusted'.
        pname : str, optional
            The name of the portfolio. The default is 'Port'.
        pcolname : str, optional
            Name of the portfolio price column. If it set to `None` then 
            `pcolname=pname`. The default is `None`.
        capital : float, optional
            Initial portfolio Capital in dollars. The default is 100000.

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
        
    def set_model(self, ww):
        """
        Set model parameters and evaluate the portfolio time-series.

        Parameters
        ----------
        ww : pandas.DataFrame
            Rebalance schedule. The following columns must be present: \n
                "Droll" : rolling date (rebalancing day) \n
                "Dfix" : fixing date (day when the close price are recorded) \n
                name of symbol 1 : weights for symbol 1 \n
                name of symbol 2 : weights for symbol 2 \n
                erc. \n
            The symbol name must match the symbol names from `mktdata`.
            
        Returns
        -------
        pandas.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        # Make sure that the weights are normalized
        self.ww = ww.copy()
        self.ww[self.symb] = self.ww[self.symb] \
                                 .apply(lambda x: x / x.sum(), axis=1)
        
        # Calculate
        self._port_calc()
        return self.port
    
    def get_nshares(self):
        """
        Returns the number of shares hold after each rolling date.

        Returns
        -------
        pandas.DataFrame
        """
        return self.nshares.astype('int').copy()
    
    def get_weights(self, fancy=False):
        """
        Returns the portfolio weights at each rebalancing period.
        
        Parameters
        ----------
        fancy : Boolean, optional
            * False: reports the weights in algebraic format.
            * True: reports the weights in percent rounded to 2 decimals.
            
        The default is False.

        Returns
        -------
        pandas.DataFrame
        """
        res = self.ww.copy()
        if not fancy:
            return res
        
        res[self.symb] = res[self.symb].round(4).abs() * 100
        
        return res
    
    def _port_calc(self):
        def _rank(x, ww=self.ww):
            for k in range(len(ww)):
                if x < ww.Droll[k]: return k
            return len(ww)
        
        mktdata = self.mktdata.pivot(columns='symbol', values=self.col_price)
        
        div = self.mktdata.pivot(columns='symbol', values=self.col_divd) \
            .groupby(_rank).sum()
        self.port = []
        self.nshares = []
        
        self.cash_invst = []
        self.cash_roll = []
        self.cash_divd = [0.]
        
        cap = self.capital
        for k, v in mktdata.groupby(_rank):
            if k == 0:  continue
 
            nsh = (self.ww[self.symb].iloc[k - 1] \
                   / mktdata.loc[self.ww.Dfix[k - 1]] * cap).round(0)
            self.nshares.append(nsh)
            self.port.append(v @ nsh)
            
            invst = nsh @ mktdata.loc[self.ww.Droll[k - 1]]
            dcap = cap - invst
            divd = div.iloc[k] @ nsh
            cap = self.port[-1][-1] + divd  + dcap
            
            self.cash_invst.append(invst)
            self.cash_roll.append(dcap)
            self.cash_divd.append(divd)
 
        self.port = pd.concat(self.port) \
            .pipe(pd.DataFrame, columns=[self.pcolname])

        self.nshares = pd.DataFrame(self.nshares,
                                    index=self.ww.Droll[:len(self.nshares)])
 
    def get_account(self, fancy=False):
        """
        Returns additional bookkeeping information regarding rebalancing 
        (e.g. residual cash due rounding number of shares, previous period 
         dividend cash accumulation, etc.)

        Parameters
        ----------
        fancy : Boolean, optional
            * `False`: the values are reported in unaltered algebraic format. 
            * `True` : the values are reported rounded.
            
        The default is `False`.

        Returns
        -------
        pandas.DataFrame
            Reports, for each rolling period identified by 'Droll': 

            * the number of shares hold for each symbol,
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
        fancy : Boolean, optional
           * `False`: returns in algebraic form.
           * `True`: returns in percent rounded to 2 decimals.
        The default is False.

        Returns
        -------
        pandas.DataFrame
            Each rolling period is indicated by its start date, Droll. 
            Included are: the fixing data, Dfix, and the 
            portfolio weights.
        """
        # local function
        def frr(x):
            p2 = x.p1.shift(-1)
            p2.iloc[-1] = self.port.iloc[-1,-1]
            return p2 / x.p1 - 1
        
        rww = self.ww.merge(self.port, left_on='Droll', right_index=True)\
                  .rename(columns={self.port.columns[0]: 'p1'})\
                  .assign(RR=frr)[['Droll', 'Dfix', 'RR'] + list(self.symb)]
 
        if not fancy:
            return rww
        
        rww['RR'] = rww['RR'].round(4) * 100
        rww[self.symb] = rww[self.symb].round(4).abs() * 100
        
        return rww
    