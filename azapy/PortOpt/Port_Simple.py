# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 14:04:34 2021

@author: mircea
"""
import pandas as pd
import numpy as np
import ta 
import numbers
import plotly.graph_objects as go

from azapy.util.drawdown import drawdown, max_drawdown

class Port_Simple:
    """
    Portfolio with constant initial weigths (no rebalance or dividend 
    accumulation).
    It is also a base class for all Port_XXX classes
    """
    def __init__(self, rprice, symb=None, sdate=None, edate=None, 
                 col='adjusted', pname='Simple', pcolname=None, 
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
        if sdate is None:
            sdate = rprice.groupby('symbol').apply(lambda x: x.index[0]).max()
        if edate is None:
            edate = rprice.groupby('symbol').apply(lambda x: x.index[-1]).min()
        if symb is None: symb = rprice.symbol.unique()
            
        self.mktdata = rprice[(rprice.index >= sdate) & 
                              (rprice.index <= edate) &
                              rprice.symbol.isin(symb)].copy()
        self.symb = pd.Series(symb)
        self.sdate = pd.to_datetime(sdate)
        self.edate = pd.to_datetime(edate)
        self.col = col
        self.pname = pname
        self.pcolname = self.pname if pcolname is None else pcolname
        self.capital = capital
        
        self.ww = None
        self.port = None
        
    def get_port(self, ww=None):
        """
        Evaluate the portfolio timeseries.

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
            self.ww = pd.Series(1., index=self.symb.size)
        elif isinstance(ww, pd.core.series.Series):
            self.ww = ww
        else:
            self.ww = pd.Series(ww, index=self.symb)
 
        assert self.ww.size == self.symb.size, \
            f"ww size must have = number of symbols {self.symb.size}"
        assert np.all(self.ww >= 0.), \
            "ww elements must be >= 0"
        wws = self.ww.sum()
        assert wws > 0, "at least one ww element must be > 0"
        self.ww = self.ww / wws
        
        self._port_calc()
        return self.port
    
    def _port_calc(self):
        mktdata = self.mktdata.pivot(columns='symbol', values=self.col)
        divid = self.ww / mktdata.iloc[0] * self.capital

        self.port = pd.DataFrame(mktdata @ divid, columns = [self.pcolname])
        
    def get_mktdata(self):
        """
        Return the actual MkT data relevant for the portfolio evaluation.
        
        Returns
        -------
        pd.DataFrame

        """
        return self.mktdata
        
    def port_view(self, emas=[30, 200], boillinger=False, 
                  view=True, fancy=False):
        """
        Plot the portfolio time series together with optional technical 
        indicators.

        Parameters
        ----------
        emas : list of int, optional
            Values for EMA durations.The default is [30, 200].
        boillinger : boolean, optional
            if set True it adds the billinger bands. The default is False.
        view : boolean, optional
            If False it suppres the ploting to terminal. The default 
            is True.
        fancy : boolena, optional
            False : it uses the pandas plot (matplotlib) capabilities.
            
            True : it uses plotly library for interactive timeseries view.
            
            The default is False.

        Returns
        -------
        df : panda.DataFrame
            The DataFrame with the timeseries included in plot.

        """
        df = self.port.copy()

        if isinstance(emas, list):
            for ema in emas:
                name = 'EMA' + str(ema)
                idx = ta.trend.EMAIndicator(df[self.pcolname], window=ema, 
                                            fillna=True)
                df[name] = idx.ema_indicator()
                
        if boillinger:
            idx = ta.volatility.BollingerBands(df[self.pcolname], fillna=True)
            df['B_upper'] = idx.bollinger_hband()
            df['B_mvag'] = idx.bollinger_mavg()
            df['B_lower'] = idx.bollinger_lband()
        
        if view: 
            if fancy: 
                self._view_plotly(df)
            else:
                df.plot()
                
        return df
    
    def port_view_all(self, sdate=None, edate=None, view=True, 
                      fancy=False):
        """
        Plot the portfolio and its component timeseries.

        Parameters
        ----------
        sdate : datetime, optional
            Start date of ploted timeseries. If it is set to None
            then the sdate is set to the erliest date in the timeseries. 
            The default is None.
        edate : TYPE, optional
            End date of ploted timesereis. If it set to None then the edate
            is set to the most recent date of the timesereis. 
            The default is None.
        view : boolean, optional
            If set to True then the plot is printed to the terminal.
            The default is True.
        fancy : boolean, optional
            False : it uses the pandas plot (matplotlib) capabilities.
            
            True : it uses plotly library for interactive timeseries view.
            
            The default is False.

        Returns
        -------
        df : pandas.DataFrame
            A Data Frame containing the timeseries.

        """
        if sdate is None: 
            sdate = self.sdate
        if edate is None: 
            edate = self.edate
        
        df = self.mktdata.pivot(columns='symbol', values=self.col) 
        df = df.iloc[(df.index >= sdate) & (df.index <= edate)] \
                .merge(self.port, on='date') \
                .apply(lambda x: x / x[0])
 
        if view: 
            if fancy: 
                self._view_plotly(df)
            else:
                df.plot()
                
        return df
        
    def port_drawdown(self, top=5, fancy=False):
        """
        Compute the portfolio drawdowns.

        Parameters
        ----------
        top : int, optional
            The number of largest drawdown that will be reported. 
            The default is 5.
        fancy : boolean, optional
            False : The drawdowns values are reported in unaltered 
            algebraic format.
            
            True : The drawdowns values are reported in percents 
            rounded to 2 decimals.
            
            The default is False.

        Returns
        -------
        res : panda.DataFrame
            Table with drawdown events with columns: \n
                "DD" : drawdown rate \n
                "Date" : recorded date of the drawdown \n
                "Star" : start date of the drawdown \n
                "End" : end date of the drawdown

        """
        res = drawdown(self.port, col=self.pcolname, top=top)
        if not fancy: return res
        
        res.DD = res.DD.round(4) * 100
        return res
    
    def port_perf(self, fancy=False):
        """
        Brief description of porfolio and its components performances
        in terms of average historical rate of returns and maximum drawdowns.

        Parameters
        ----------
        fancy : boolean, optional
            False : The rate of returns and drawdown values are reported 
            in unaltered algebraic format.
            
            True : The rate of returns and drawdown values are reported in 
            percents rounded to 2 decimals.
            
            The default is False.

        Returns
        -------
        pandas.DataFrame
            Table of porfolio and its components performance with the 
            following columns: \n
                "RR " : rate of returns \n
                "DD" : maximum rate of drawdown \n
                "DD_date" : recorder date of maximum drawdown \n
                "DD_start" : start date of maximum drawdown \n
                "DD_end" : end date of maximum drawdown
    
        """
        # local function
        def rinfo(df, col):
            rr = (df[col][-1] / df[col][0]) ** (254. / len(df)) - 1
            dv, dd, ds, de = max_drawdown(df, col=col)
            
            return pd.DataFrame([[rr, dv, dd, ds, de]],
                        columns=['RR', 'DD', 'DD_date', 'DD_start', 'DD_end' ])
        
        res = rinfo(self.port, self.pcolname)
        res.index = pd.Index([self.pname], name='symbol')

        res2 = self.mktdata.groupby('symbol').apply(rinfo, 
                                                    col=self.col).droplevel(1)
 
        res = pd.concat([res, res2])
        
        if not fancy: 
            return res
        
        res.RR = res.RR.round(4) * 100
        res.DD = res.DD.round(4) * 100
        return res
        
    def port_annual_returns(self, fancy=False):
        """
        Porfolio historical calendar annual rate of returns

        Parameters
        ----------
        fancy : boolean, optional
            False : The rates are reported in unaltered 
            algebraic format.
            
            True :The rates are reported in percents rounded to 2 decimals 
            and presented is color style.
            
            The default is False.

        Returns
        -------
        pandas.DataFrame

        """
        # local function
        def rrate(df):
            return df[self.pcolname][-1] / df[self.pcolname][0] - 1
        
        res = self.port.resample('A', convention='end').apply(rrate)
        res.index = res.index.year
        res.index.name = 'year'
        res.rename('RR', inplace=True)
        
        if not fancy: 
            return res
            
        return pd.DataFrame(res.round(4)) \
                 .style.format("{:.2%}") \
                 .applymap(self._color_negative_red)
        
    
    def port_monthly_returns(self, fancy=False):
        """calendaristic 
        Portfolio historical calendar monthly rate of returns.

        Parameters
        ----------
        fancy : boolean, optional
            False : The rates are reported in unaltered 
            algebraic format.
            
            True :The rates are reported in percents rounded to 2 decimals 
            and presented is color style.
            
            The default is False.

        Returns
        -------
        pandas.DataFrame

        """
        def rrate(df):
            return df[self.pcolname][-1] / df[self.pcolname][0] - 1
        
        res = self.port.resample('M', convention='end').apply(rrate)
        
        if not fancy: 
            return res
        
        res.index = pd.MultiIndex.from_arrays([res.index.year.to_list(), 
                                               res.index.month.to_list()],
                                              names= ['year', 'month'])
    
        res = res.reset_index(level=0).pivot(columns='year', 
                                             values=0).round(4) 
 
        return res.style.format("{:.2%}").applymap(self._color_negative_red)
 
    
    def _color_negative_red(self, val):
        if isinstance(val, numbers.Number) and (val < 0):
            return 'color: red'
        else:
            return 'color: black'
        
    def _view_plotly(self, df):
        # Create figure
        fig = go.Figure()

        for col in df.columns:
            fig.add_trace(go.Scatter(x=list(df.index), y=list(df[col]), 
                                     name=col, line=dict(width=1.5)))

        # Set title
        fig.update_layout(
            title_text="Port performance"
        )

        # Add range slider
        fig.update_layout(
            xaxis=dict(
                rangeselector=dict(
                    buttons=list([
                        dict(count=1, label="1m", 
                             step="month", stepmode="backward"),
                        dict(count=6, label="6m", 
                             step="month", stepmode="backward"),
                        dict(count=1, label="YTD", 
                             step="year", stepmode="todate"),
                        dict(count=1, label="1y", 
                             step="year", stepmode="backward"),
                        dict(step="all")
                    ])
                ),
                rangeslider=dict(
                    visible=True,
                    thickness = 0.05
                ),
                autorange = False,
                range=[df.index[0], df.index[-1]],
                fixedrange= False,
                type="date"
            )
        )

        fig.update_layout(
            yaxis=dict(
                autorange = False,
                range=[df.min().min() * 0.95, df.max().max() * 1.05],
                fixedrange= False
            )
        )

        fig.update_layout(
            autosize=False,
            width=900,
            height=550,
            margin=dict(l=50, r=50, b=25, t=100, pad=1),
            paper_bgcolor="LightSteelBlue",
        )

        fig.show()