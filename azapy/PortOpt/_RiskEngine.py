import numpy as np
import pandas as pd

from azapy.MkT.MkTcalendar import NYSEgen

class _RiskEngine():
    """
    Base class. Derive class needs to implement\n
        getWeights \n

    Methods:
        * getPositions
        * set_rrate
        * set_mktdata
    """

    def __init__(self, mktdata=None, colname='adjusted',
                 freq='Q', hlength=3.25, calendar=None):
        """
        Constructor

        Parameters
        ----------
        `mktdata` : pandas.DataFrame, optional
            Historic daily market data for portfolio components in the format
            returned by azapy.mktData function. The default is `None`.
        `colname` : string, optional
            Name of the price column from mktdata used in the weights
            calibration. The default is 'adjusted'.
        `freq` : str, optional
            Rate of returns horizon in number of business day. it could be
            'Q' for quarter or 'M' for month. The default is 'Q'.
        `hlength` : float, optional
            History length in number of years used for calibration. A
            fractional number will be rounded to an integer number of months.
            The default is 3.25
        `calendar` : numpy.busdaycalendar, optional
            Business days calendar. If is it None then the calendar will be set
            to NYSE business calendar.
            The default is `None`.

        Returns
        -------
        The object.
        """
        if mktdata is not None:
            self.set_mktdata(mktdata, colname, freq, hlength)

    def set_rrate(self, rrate):
        """
        Sets portfolio components historical rates of returns in the format
        "date", "symbol1", "symbol2", etc.

        Parameters
        ----------
        `rrate` : pandas.DataFrame
            The portfolio components historical rates of returns.
            If it is not None, it will overwrite the rrate computed in the
            constructor from mktdata. The default is None.
        Returns
        -------
        None.
        """
        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean()
        self.rrate = rrate

    def set_mktdata(self, mktdata, colname='adjusted',
                    freq='Q', hlength=3.25, calendar=None):
        """
        Sets historical market data. It will overwrite the choice made in the
        constructor.

        Parameters
        ----------
        `mktdata` : pandas.DataFrame
            Historic daily market data for portfolio components in the format
            returned by azapy.mktData function.
        `colname` : str, optional
            Name of the price column from mktdata used in the weights
            calibration. The default is 'adjusted'.
        `freq` : string, optional
            Rate of returns horizon in number of business day. it could be
            'Q' for quarter or 'M' for month. The default is 'Q'.
        `hlength` : float, optional
            History length in number of years used for calibration. A
            fractional number will be rounded to an integer number of months.
            The default is 3.25
        `calendar` : numpy.busdaycalendar, optional
            Business days calendar. If is it None then the calendar will be set
            to NYSE business calendar.
            The default is `None`.

        Returns
        -------
        None.
        """
        if mktdata is None:
            return

        self.mktdata = mktdata

        periods = 63 if freq == 'Q' else 21

        if calendar is None:
            calendar = NYSEgen()

        edate = mktdata.index[-1]
        sdate = edate - pd.offsets.DateOffset(months=round(hlength * 12, 0))
        sdate = np.busday_offset(sdate.date(),
                                 0, roll='backward',  busdaycal=calendar)
        rrate = mktdata.pivot(columns='symbol', values=colname)[sdate:]
        self.last_data = rrate.iloc[-1]
        rrate = rrate.pct_change(periods=periods).dropna()

        self.set_rrate(rrate)

    def getPositions(self, nshares=None, cash=0., ww=None):
        """
        Computes the number of shares according to the weights

        Parameters
        ----------
        `nshares` : pandas.Series, optional
            Number of shares per portfolio component. A missing component
            entry will be considered 0. A `None` value assumes that all
            components entries are 0. The name of the components must be
            present in the mrkdata. The default is `None`.
        `cash` : float, optional
            Additional cash to be considered in the overall capital. A
            negative entry assumes a reduction in the total capital
            available for rebalance. The default is 0.
        `ww` : pandas.Series, optional
            External portfolio weights. If it not set to `None` these
            weights will overwrite the calibrated weights.
            The default is `None`.

        Returns
        -------
        pandas.DataFrame:
        the rolling information.

        Columns:

            - 'old_nsh' :
                the initial number of shares per portfolio
                component as well as additional cash position. These are
                present in the input.
            - 'new_nsh' :
                the new number of shares per component plus the
                residual cash (due to the rounding to an integer number of
                shares). A negative entry means that the investor needs to
                add more cash to cover for the number of share
                roundups. It has a small value.
            - 'diff_nsh' :
                the number of shares that needs to be
                both/sold in order to rebalance the portfolio positions.
            - 'weights' :
                portfolio weights used for rebalanceing. The 'cash'
                entry is the new portfolio value.
            - 'prices' :
                the share prices used for rebalances evaluations.

        Note: Since the prices are closing prices, the rebalance can be
        executed next business. Additional cash slippage may occur due
        to share price differential between the previous day closing and
        execution time.

        """
        ns = pd.Series(0, index=self.rrate.columns)
        if nshares is not None:
            ns = ns.add(nshares, fill_value=0)

        if len(self.rrate.columns) != len(ns):
            raise ValueError("nshares must by a subset "
                             f"of {self.rrate.columns}")

        pp = self.last_data[self.rrate.columns.to_list()]

        cap = ns.dot(pp) + cash
        if ww is None:
            ww = self.getWeights()
        else:
            ww0 = pd.Series(0, index=self.rrate.columns)
            ww = ww0.add(ww, fill_value=0)
            if len(self.rrate.columns) != len(ns):
                raise ValueError("ww: wrong component names - "
                                 f"must be among {self.rrate.columns}")

        newns = (ww / pp * cap).round(0)
        newcash = cap - newns.dot(pp)

        ns['cash'] = cash
        newns['cash'] = newcash
        ww['cash'] = cap - newcash
        pp['cash'] = 1.

        res = pd.DataFrame({'old_nsh': ns,
                            'new_nsh': newns,
                            'diff_nsh': newns - ns,
                            'weights' : ww.round(6),
                            'prices': pp})

        return res

    # need to be implemented by the derive classes
    def getWeights(self):
        pass
