import pandas as pd
import numpy as np
import ta
import numbers
import plotly.graph_objects as go

from azapy.util.drawdown import drawdown, max_drawdown

class Port_Simple:
    """
    Back testing the Buy and Hold portfolio.

    Methods:
        * set_model
        * get_port
        * get_mktdata
        * port_view
        * port_view_all
        * port_drawdown
        * port_perf
        * port_annual_returns
        * port_monthly_returns
    """

    def __init__(self, mktdata, symb=None, sdate=None, edate=None,
                 col='adjusted', pname='Port', pcolname=None,
                 capital=100000):
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
            should be included in `mktdata`. If set to `None` the `symb` will 
            be  set to the full set of symbols included in `mktdata`. 
            The default is `None`.
        sdate : date like, optional
            Start date for historical data. If set to `None` the `sdate` will
            be set to the earliest date in `mktdata`. The default is `None`.
        edate : date like, optional
            End date for historical dates and so the simulation. Must be
            greater than `sdate`. If it is `None` then `edate` will be set
            to the latest date in `mktdata`. The default is `None`.
        col : str, optional
            Name of column in the mktdata DataFrame that will be considered
            for portfolio aggregation. The default is 'adjusted'.
        pname : str, optional
            The name of the portfolio. The default is 'Port'.
        pcolname : str, optional
            Name of the portfolio price column. If it set to `None` that
            `pcolname=pname`. The default is `None`.
        capital : float, optional
            Initial portfolio Capital in dollars. The default is 100000.

        Returns
        -------
        The object.
        """
        if isinstance(mktdata, list):
            mktdata = pd.concat(mktdata, axis=1, join='inner') \
                .melt(var_name='symbol', value_name=col, ignore_index=False)

        # set sdate
        mkt_sdate = mktdata.groupby('symbol') \
                           .apply(lambda x: x.index[0]).max()
        if sdate is None:
            self.sdate = mkt_sdate
        else:
            v_date = pd.to_datetime(sdate)
            if v_date < mkt_sdate:
                self.sdate = mkt_sdate
            else:
                self.sdate = v_date
        # set edate
        mkt_edate = mktdata.groupby('symbol') \
                           .apply(lambda x: x.index[-1]).min()
        if edate is None:
            self.edate = mkt_edate
        else:
            v_date = pd.to_datetime(edate)
            if v_date > mkt_edate:
                self.edate = mkt_edate
            else:
                self.edate = v_date
        # set symb
        mkt_symb = mktdata.symbol.unique()
        if symb is None:
            self.symb = mkt_symb
        else:
            if not all(sy in mkt_symb for sy in symb):
                raise ValueError(f"symb list {symb} incompatible"
                                 f" with mktdata {mkt_symb}")

            self.symb = pd.Series(symb)
        # set mktdata
        self.mktdata = mktdata[(mktdata.index >= self.sdate) &
                               (mktdata.index <= self.edate) &
                               mktdata.symbol.isin(self.symb)].copy()

        self.col = col
        self.pname = pname
        self.pcolname = self.pname if pcolname is None else pcolname
        self.capital = capital

        self.ww = None
        self.port = None

    def set_model(self, ww=None):
        """
        Set model parameters and evaluate the portfolio time-series.

        Parameters
        ----------
        ww : list (also numpy.array or pandas.Series), optional
            List of weights. If it is pandas.Series the index should match
            the basket `symb`. Otherwise the weights are considered in the 
            `symb` order. If it is set to `None` than `ww` will be set 
            to equal weights.
            The default is `None`.

        Returns
        -------
        pandas.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        # Validate the weights
        if ww is None:
            self.ww = pd.Series(1., index=self.symb)
        elif isinstance(ww, pd.core.series.Series):
            self.ww = ww
        else:
            self.ww = pd.Series(ww, index=self.symb)

        if self.ww.size != self.symb.size:
            raise ValueError(f"wrong ww size - must be {self.symb.size}")
        if np.any(self.ww < 0.):
            raise ValueError("ww elements must be >= 0")
        wws = self.ww.sum()
        if wws <= 0.:
            raise ValueError("At least one ww element must be > 0")
        self.ww = self.ww / wws

        # Calculate
        self._port_calc()
        return self.port

    def _port_calc(self):
        mktdata = self.mktdata.pivot(columns='symbol', values=self.col)
        divid = self.ww / mktdata.iloc[0] * self.capital

        self.port = pd.DataFrame(mktdata @ divid, columns = [self.pcolname])

    def get_port(self):
        """
        Returns the portfolio time-series.

        Returns
        -------
        pandas.DataFrame
        """
        return self.port.copy()

    def get_mktdata(self):
        """
        Returns the actual market data used for portfolio evaluations.

        Returns
        -------
        pandas.DataFrame
        """
        return self.mktdata.copy()

    def port_view(self, emas=[30, 200], bollinger=False, fancy=False,
                  saveto=None):
        """
        Plot the portfolio time series together with optional technical
        indicators.

        Parameters
        ----------
        emas : list of int, optional
            List of EMA durations. The default is [30, 200].
        bollinger : Boolean, optional
            If set `True` it adds the Bollinger bands. The default is `False`.
        fancy : Boolean, optional
            `False` : it uses the matplotlib capabilities.

            `True` : it uses plotly library for interactive time-series view.

            The default is `False`.
        saveto : str, optional
            The name of the file where to save the plot. The default is `None`.

        Returns
        -------
        pandas.DataFrame
            Contains the time-series included in plot.
        """
        df = self.port.copy()

        if isinstance(emas, list):
            for ema in emas:
                name = 'EMA' + str(ema)
                idx = ta.trend.EMAIndicator(df[self.pcolname], window=ema,
                                            fillna=True)
                df[name] = idx.ema_indicator()

        if bollinger:
            idx = ta.volatility.BollingerBands(df[self.pcolname], fillna=True)
            df['B_upper'] = idx.bollinger_hband()
            df['B_mvag'] = idx.bollinger_mavg()
            df['B_lower'] = idx.bollinger_lband()

        if fancy:
            fig = self._view_plotly(df)
            if saveto is not None:
                fig.write_image(saveto)
        else:
            if saveto is None:
                df.plot()
            else:
                df.plot().get_figure().savefig(saveto)

        return df

    def port_view_all(self, sdate=None, edate=None, componly=False,
                      fancy=False, saveto=None):
        """
        Plot the portfolio and its component time-series in a relative bases.

        Parameters
        ----------
        sdate : date like, optional
            Start date of plotted time-series. If it is set to `None`
            then the `sdate` is set to the earliest date in the time-series.
            The default is `None`.
        edate : date like, optional
            End date of plotted time-series. If it set to `None` then 
            the `edate
            is set to the most recent date of the time-series.
            The default is `None`.
        componly : Boolean, optional
            `True` : only the portfolio components time-series are plotted.

            `False`: the portfolio and its components times-series are plotted.

            The default is `True`.
        fancy : Boolean, optional
            `False` : it uses the pandas plot (matplotlib) capabilities.

            `True` : it uses plotly library for interactive time-series view.

            The default is `False`.
        saveto : str, optional
            The name of the file where to save the plot. The default is `None`.

        Returns
        -------
        pandas.DataFrame
            A Data Frame containing the time-series.
        """
        if sdate is None:
            sdate = self.sdate
        if edate is None:
            edate = self.edate
            
        sdate = pd.to_datetime(sdate)
        edate = pd.to_datetime(edate)

        df = self.mktdata.pivot(columns='symbol', values=self.col)
        df = df.iloc[(df.index >= sdate) & (df.index <= edate)]
        if not componly:
            df = df.merge(self.port, on='date')
        df = df.apply(lambda x: x / x[0])

        if fancy:
            fig = self._view_plotly(df)
            if saveto is not None:
                fig.write_image(saveto)
        else:
            if saveto is None:
                df.plot()
            else:
                df.plot().get_figure().savefig(saveto)

        return df

    def port_drawdown(self, top=5, fancy=False):
        """
        Computes the portfolio drawdowns.

        Parameters
        ----------
        top : int, optional
            The number of largest drawdowns that will be reported.
            The default is 5.
        fancy : Boolean, optional
            `False` : The drawdowns values are reported in unaltered
            algebraic format.

            `True` : The drawdowns values are reported in percent
            rounded to 2 decimals.

            The default is `False`.

        Returns
        -------
        res : panda.DataFrame
            Table of drawdown events. Columns: \n
                'DD' : drawdown rate \n
                'Date' : recorded date of the drawdown \n
                'Star' : start date of the drawdown \n
                'End' : end date of the drawdown
        """
        res = drawdown(self.port, col=self.pcolname, top=top)
        if not fancy: return res

        res.DD = res.DD.round(4) * 100
        return res

    def port_perf(self, componly=False, fancy=False):
        """
        Brief description of portfolio and its components performances
        in terms of average historical rate of returns and maximum drawdowns.

        Parameters
        ----------
        componly : Boolean, optional
            If `True`, only the portfolio components maximum drawdowns
            are reported. The default is `False`.
        fancy : Boolean, optional
            * `False` : The rate of returns and drawdown values are reported
            in unaltered algebraic format.

            * `True` : The rate of returns and drawdown values are reported
            in  percent rounded to 2 decimals.

            The default is `False`.

        Returns
        -------
        pandas.DataFrame
            Performance information.
            Columns:

            * 'RR' : rate of returns
            * 'DD' : maximum rate of drawdown
            * 'Beta' : abs(RR/DD)
            * 'DD_date' : recorder date of maximum drawdown
            * 'DD_start' : start date of maximum drawdown
            * 'DD_end' : end date of maximum drawdown
        """
        # local function
        def rinfo(df, col):
            rr = (df[col][-1] / df[col][0]) ** (252. / len(df)) - 1
            dv, dd, ds, de = max_drawdown(df, col=col)

            return pd.DataFrame([[rr, dv, np.abs(rr / dv), dd, ds, de]],
                columns=['RR', 'DD', 'Beta', 'DD_date', 'DD_start', 'DD_end' ])

        res = self.mktdata.groupby('symbol') \
            .apply(rinfo, col=self.col) \
            .droplevel(1) \
            .sort_values(by='Beta', ascending=False)

        if not componly:
            res2 = rinfo(self.port, self.pcolname)
            res2.index = pd.Index([self.pname], name='symbol')
            res = pd.concat([res2, res])

        if not fancy:
            return res

        res.RR = res.RR.round(4) * 100
        res.DD = res.DD.round(4) * 100
        return res

    def port_annual_returns(self, withcomp=False, componly=False, fancy=False):
        """
        Portfolio annual (calendar) rates of returns.

        Parameters
        ----------
        withcomp : Boolean, optional
            If `True`, adds the portfolio components annual returns to the
            report. The default is `False`.
        componly : Boolean, optional
            If `True`, only the portfolio components annual returns are 
            reported.
            The flag is active only if `withcomp=True`. The default is `False`.
        fancy : Boolean, optional
            `False` : The rates are reported in unaltered
            algebraic format.

            `True` :The rates are reported in percent rounded to 2 decimals
            and presented is color style.

            The default is `False`.

        Returns
        -------
        pandas.DataFrame
        """
        # local function
        def frrate(df):
            return df[-1] / df[0] - 1

        if withcomp:
            zz = self.mktdata.pivot(columns='symbol', values=self.col) \
                .loc[self.port.index]
            if not componly:
                zz = self.port.merge(zz, on='date')
            res = zz.resample('A', convention='end').apply(frrate)
        else:
            res = self.port.resample('A', convention='end').apply(frrate)
        res.index = res.index.year
        res.index.name = 'year'

        if not fancy:
            return res

        return pd.DataFrame(res.round(4)) \
                 .style.format("{:.2%}") \
                 .applymap(self._color_negative_red)


    def port_monthly_returns(self, withcomp=False, componly=False,
                             fancy=False):
        """
        Portfolio monthly (calendar) rate of returns.

        Parameters
        ----------
        withcomp : Boolean, optional
            If `True`, adds the portfolio components monthly returns to the
            report. The default is `False`.
        componly : Boolean, optional
            If `True`, only the portfolio components monthly returns are
            reported. The flag is active only if `withcomp=True`.
            The default is `False`.
        fancy : Boolean, optional
            `False` : The rates are reported in unaltered
            algebraic format.

            `True` :The rates are reported in percent rounded to 2 decimals
            and presented is color style.

            The default is `False`.

        Returns
        -------
        pandas.DataFrame
        """
        def frrate(df):
            return df[-1] / df[0] - 1

        # res = self.port.resample('M', convention='end').apply(frrate)
        if withcomp:
            zz = self.mktdata.pivot(columns='symbol', values=self.col) \
                .loc[self.port.index]
            if not componly:
                zz = self.port.merge(zz, on='date')
            res = zz.resample('M', convention='end').apply(frrate)
        else:
            res = self.port.resample('M', convention='end').apply(frrate)

        if not fancy:
            return res

        res.index = pd.MultiIndex.from_arrays([res.index.year.to_list(),
                                               res.index.month.to_list()],
                                              names= ['year', 'month'])

        if not withcomp:
            res = res.reset_index(level=0).pivot(columns='year',
                                                 values=self.pcolname).round(4)

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

        return fig
