# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 15:06:34 2021

@author: mircea
"""

import pandas as pd

from .Port_Simple import Port_Simple

class Port_Weighted(Port_Simple):
    """
    Portfolio with custom scheduled weights.
    Iherited from azapy.Port_simple \n
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
    def __init__(self, rprice, symb=None, sdate=None, edate=None, 
                 col='close', pname='Weighted', pcolname=None,
                 capital=100000):
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
            for portfolio agregation.The default is 'adjusted'.
        pname : string, optional
            The name of the portfolio. The default is 'Simple'.
        pcolname : string, optional
            Name of the portfolio price column. If it set to None that 
            pcolname=pname. The default is None.
        capital : float, optional
            Initial portfolio Capital in dollars. The default is 100000.

        Returns
        -------
        The object.

        """
        super().__init__(rprice=rprice, symb=symb, 
                         sdate=sdate, edate=edate, 
                         col=col, pname=pname, 
                         pcolname=pcolname, capital=capital)
        self.nshares = None
        
        self.cash_invst = None
        self.cash_roll = None
        self.cash_divd = None
        
    def get_port(self, ww):
        """
        Evaluate portfolio timeseries.

        Parameters
        ----------
        ww : pd.DataFrame
            Rebalace schedule. The following columns must be present: \n
                "Droll" : rolling date (rebalancing day) \n
                "Dfix" : fixing date (day when the close price are recorded) \n
                name of symbol 1 : weights for symbol 1 \n
                name of symbol 2 : weights for symbol 2 \n
                erc. \n
            The symbol name must match the symbol names from rprice.
            
        Returns
        -------
        pd.DataFrame
            The portfolio timeseries in the format "date", "pcolname".

        """
        self.ww = ww.copy()
        self.ww[self.symb] = self.ww[self.symb] \
            .apply(lambda x: x/x.sum(), axis=1)
        
        self._port_calc()
        return self.port
    
    def get_nshares(self):
        """
        Number of shares hold after each rebalance.

        Returns
        -------
        pandas.DataFrame

        """
        return self.nshares.astype('int')
    
    def get_weights(self):
        """
        Returns the portfolio weights at each rebalancing period

        Returns
        -------
        pandas.DataFrame

        """
        
        return self.ww
    
    def _port_calc(self):
        def _rank(x, ww=self.ww):
            for k in range(len(ww)):
                if x <= ww.Droll[k]: return k
            return len(ww)
        
        mktdata = self.mktdata.pivot(columns='symbol', values=self.col)
        
        div = self.mktdata.pivot(columns='symbol', values='divd') \
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
                                    index=self.ww.Droll[:-1])
 
    def get_account(self, fancy=False):
        """
        Returns additional bookkeeping information regarding rebalancing 
        (e.g. residual cash due rounding number of shares, previouse period 
         dividend cash accumulation, etc.)

        Parameters
        ----------
        fancy : boolean, optional
            False: returns a simple pandas.DataFrame
            
            True : returns a pandas.DataFrame with rounded values.
            
           The default is False.

        Returns
        -------
        pandas.DataFrame
            For each rolling period identified by "Droll" (rolling date) index,
            reports: \n
                for each symbol : the number of shares hold \n
                "cash_invst" : cash invested at the begining of 
                period \n
                "cash_roll" : cash rolled to the next period \n
                "cash_divd" : cash dividend accumulated in the 
                previouse period \n
            Note: The capital at the begining of the period is 
            cash_invst + cash_roll. It is also equal to the previous period: 
            value of the shares on the fixing date + cash_roll + cash_divd.
            There are 2 sources for the cash_roll. The roundup to integer 
            number of shares and the shares close price differences between 
            the fixing (computation) and rolling (execution) dates. It could
            be positive or negative. The finance of the cash_roll during 
            each rolling period is assumed  to be done separately by the 
            investor.
        """
        acc_tab = self.nshares.copy()
        acc_tab['cash_invst'] = self.cash_invst
        acc_tab['cash_roll'] = self.cash_roll
        acc_tab['cash_divd'] = self.cash_divd[:-1]
        
        if not fancy: return acc_tab
        
        acc_tab = acc_tab.round(2)
        acc_tab[self.symb] = acc_tab[self.symb].astype('int')
        
        return acc_tab
