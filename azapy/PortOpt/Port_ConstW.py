import pandas as pd
import numpy as np

from .Port_Rebalanced import Port_Rebalanced
from azapy.MkT.MkTcalendar import NYSEgen
from azapy.util.schedule import schedule_roll

class Port_ConstW(Port_Rebalanced):
    """
    Backtesting portfolio with constant weights, periodically rebalanced.

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
                 pname='Port', pcolname=None, capital=100000,
                 schedule=None,
                 freq='Q', noffset=-3, fixoffset=-1, calendar=None):
        """
        Constructor

        Parameters
        ----------
        `mktdata` : `pandas.DataFrame`
            MkT data in the format "symbol", "date", "open", "high", "low",
            "close", "volume", "adjusted", "divd", "split" (e.g. as returned
            by `azapy.readMkT` function).
        `symb` : list, optional
            List of symbols for the basket components. All symbols MkT data
            should be included in mktdata. If set to `None` the `symb` will be
            set to include all the symbols from `mktdata`. The default
            is `None`.
        `sdate` : date like, optional
            Start date for historical data. If set to `None` the `sdate` will
            be set to the earliest date in mktdata. The default is `None`.
        `edate` : date like, optional
            End date for historical dates and so the simulation. Must be
            greater than  `sdate`. If it is `None` then `edat`e will be set
            to the latest date in mktdata. The default is `None`.
        `col_price` : str, optional
            Column name in the mktdata DataFrame that will be considered
            for portfolio aggregation. The default is `'close'`.
        `col_divd` :  str, optional
            Column name in the mktdata DataFrame that holds the dividend
            information. The default is `'dvid'`.
        `col_ref` : str, optional
            Column name in the mktdata DataFrame that will be used as a price
            reference for portfolio components. The default is `'adjusted'`.
        `pname` : str, optional
            The name of the portfolio. The default is `'Port'`.
        `pcolname` : str, optional
            Name of the portfolio price column. If it set to `None` that
            `pcolname=pname`. The default is `None`.
        `capital` : float, optional
            Initial portfolio Capital in dollars. The default is `100000`.
        `schedule` : `pandas.DataFrame`, optional
            Rebalancing schedule, with columns for 'Droll' rolling date and
            'Dfix' fixing date. If it is `None` than the schedule will be set
            using the `freq`, `noffset`, `fixoffset` and `calendar`
            information. The default is `None`.
        `freq` : str, optional
            rebalancing frequency. It can be 'Q' for quarterly or 'M' for
            monthly rebalancing, respectively. It is relevant only is schedule
            is `None`. The default is `'Q'`.
        `noffset` : int, optional
            Number of business days offset for rebalancing date 'Droll'
            relative to the end of the period (quart or month). A positive
            value add business days beyond the calendar end of the period while
            a negative value subtracts business days. It is relevant only is
            schedule is `None`. The default is `-3`.
        `fixoffset` : int, optional
            Number of business day offset of fixing date 'Dfix' relative to
            the rebalancing date 'Droll'. It can be 0 or negative. It is
            relevant only is schedule is `None`. The default is `-1`.
        `calendar` : `numpy.busdaycalendar`, optional
            Business calendar. If it is `None` then it will be set to NYSE
            business calendar. The default
            vale is `None`.

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
        ww : list (also `numpy.array` or `pandas.Series`), optional

            List of weights. If it is `pandas.Series` the index should match
            the basket `symb`. Otherwise, the weights are considered in the 
            `symb` order. 
            If it is set to `None` than `ww` will be set to equal weights.
            The default is `None`.

        Returns
        -------
        `pands.DataFrame`
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

        if _ww.size != self.symb.size:
            raise ValueError(f"ww wrong size, it must be {self.symb.size}")

        if np.any(_ww < 0.):
            raise ValueError("All ww elements must be >= 0")

        wws = _ww.sum()
        if wws <= 0.:
            raise ValueError("At least one ww element must be > 0")

        self.ww = self.schedule
        for sy in self.symb:
            self.ww[sy] = _ww[sy] / wws
