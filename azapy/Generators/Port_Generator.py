import numpy as np
import pandas as pd

from .Port_ConstW import Port_ConstW

class Port_Generator(Port_ConstW):
    """
    Backtesting portfolio with weights genrated by a model
    
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
        `col_calib` : `str`, optional;
            Column name used for historical weights calibrations. 
            The default is `'adjusted'`.
        `pname` : `str`, optional;
            The name of the portfolio. The default is `'Port'`.
        `pcolname` : `str`, optional;
            Name of the portfolio price column. If it set to `None` that
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
        
    def set_model(self, wwModel, rtype='Sharpe',
                  mu=None, mu0=0, aversion=None, ww0=None, 
                  verbose=False):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        `rtype` : `str`, optional;
            Optimization type. Possible values: \n
                `'Risk'` : optimal-risk portfolio for targeted expected rate of 
                return.\n
                `'Sharpe'` : Sharpe-optimal portfolio - maximization solution.\n
                `'Sharpe2'` : Sharpe-optimal portfolio - minimization solution.\n
                `'MinRisk'` : minimum risk portfolio.\n
                `'RiskAverse'` : optimal-risk portfolio for a fixed 
                risk-aversion factor.\n
                `'InvNrisk'` : optimal-risk portfolio with the same risk value 
                as a benchmark portfolio (e.g., same as equal weighted 
                portfolio).\n
                `'Diverse'` : optimal-diversified portfolio for targeted
                expected rate of return (maximum of inverse 1-D).\n
                `'Diverse2'` : optimal-diversified portfolio for targeted
                expected rate of return (minmum of 1-D).\n
                `'MaxDiverse'` : maximum diversified portfolio.\n
                `'InvNdiverse'` : optimal-diversified portfolio with the same
                diversification factor as a benchmark portfolio 
                (e.g., same as equal weighted portfolio).\n
                `'InvNdrr'` : optima- diversified portfolio with the same 
                expected rate of return as a benchmark portfolio
                (e.g., same as equal weighted portfolio).\n
            The defauls is `'Sharpe'`.
        `mu` : `float`, optional;
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype='Risk'`
            The default is `None`.
        `mu0` : `float`, optional;
            Risk-free rate accessible to the investor.
            Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        `aversion` : `float`, optional;
            The value of the risk-aversion factor.
            Must be positive. Relevant only if `rtype='RiskAvers'`.
            The default is `None`.
        `ww0` : `list` (also `numpy.array` or `pandas.Series`), optional;
            Targeted portfolio weights. 
            Relevant only if `rype='InvNrisk'`.
            Its length must be equal to the number of
            symbols in rrate (mktdata). 
            All weights must be >= 0 with sum > 0.
            If it is a list or a `numpy.array` then the weights are assumed to
            by in order of `rrate.columns`. If it is a `pandas.Series` then 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (same symbols, not necessary in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.
        `hlength` : `float`, optional;
            The length in year of the historical calibration period relative
            to `'Dfix'`. A fractional number will be rounded to an integer 
            number of months. The default is `3.25` years.
        `method` : `str`, optional;
            Linear programming numerical method.
            Could be: `'ecos'`, `'highs-ds'`, `'highs-ipm'`, `'highs'`,
            `'interior-point'`, `'glpk'` and `'cvxopt'`.
            The default is `'ecos'`.
        `verbose` : Boolean, optiona;
            If it set to `True` then it will print messages when the optimal
            portfolio degenerates to a single asset portfolio as a limited 
            case. 
            The default is `False`.

         Returns
        -------
        `pandas.DataFrame`;
            The portfolio time-series in the format 'date', 'pcolname'.
        """
        self._set_alpha(alpha, coef)
        self._set_rtype(rtype)
        self.mu = mu
        self.mu0 = mu0
        self.aversion = aversion
        self.ww0 = ww0
        self.hlength = hlength
        self._set_method(method)
        self.verbose = verbose

        self._set_schedule()
        self._set_weights()
        self._port_calc()
        return self.port
    
    def set_model_X(self, wwModel):
        """
        Set model parameters and evaluate the portfolio time-series.
        
        Parameters
        ----------
        `hlength` : `float`, optional;
            The length in year of the historical calibration period relative 
            to `'Dfix'`. A fractional number will be rounded to an integer number 
            of months. The default is `3.25` years. 

        Returns
        -------
        `pandas.DataFrame`;
            The portfolio time-series in the format "date", "pcolname".
        """
        #self.hlength = hlength
        
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