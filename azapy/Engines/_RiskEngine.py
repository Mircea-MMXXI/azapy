import numpy as np
import pandas as pd

class _RiskEngine():
    """
    Base class.
    
    **Attributes**
        * status : `int` - computation status (`0` - success, any other 
          value indicates an error)
        * ww : `pandas.Series` - portfolio weights
        * name : `str` - portfolio name
    
    Note: Derive class needs to implement `getWeights` method
    
    **Methods**
        * getPositions
        * set_rrate
        * set_mktdata
    """
    def __init__(self, mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, name=None):
        """
        Constructor

        Parameters
        ----------
        mktdata : `pandas.DataFrame`, optional
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        colname : `str`, optional
            Name of the price column from mktdata used in the weights
            calibration. The default is `'adjusted'`.
        freq : `str`, optional
            Rate of return horizon. It could be 
            `'Q'` for a quarter or `'M'` for a month. The default is `'Q'`.
        hlength : `float`, optional
            History length in number of years used for calibration. A
            fractional number will be rounded to an integer number of months.
            The default is `3.25` years.
        name : `str`, optional
            Portfolio name. Deafult value is `None`

        Returns
        -------
        The object.
        """
        self._ptype_ = 'Optimizer'
        self.ww = None
        self.status = None
        
        self.colname = None
        self.freq = None
        self.hlength = None
        self.name = name
        
        self.verbose = False
        
        self.time_level1 = None
        self.time_level2 = None
        self.time_level3 = None

        self.set_mktdata(mktdata, colname, freq, hlength)


    def set_rrate(self, rrate):
        """
        Sets portfolio components historical rates of return in the format
        "date", "symbol1", "symbol2", etc.

        Parameters
        ----------
        rrate : `pandas.DataFrame`
            The portfolio components historical rates of return.
            If it is not `None`, it will overwrite the `rrate` computed in the
            constructor from `mktdata`. The default is `None`.

        Returns
        -------
        `None`
        """
        if rrate is None:
            # noting to do
            return

        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean(numeric_only=True)
        self.rrate = rrate


    def set_mktdata(self, mktdata, colname='adjusted', freq=None, hlength=None,
                    pclose=False):
        """
        Sets historical market data. It will overwrite the choice made in the
        constructor.

        Parameters
        ----------
        mktdata : `pandas.DataFrame`
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function.
        colname : `str`, optional
            Name of the price column from mktdata used in the weight's
            calibration. The default is 'adjusted'.
        freq : `str`, optional
            Rate of return horizon. It could be 
            `'Q'` for a quarter or `'M'` for a month. The default is `'Q'`.
        hlength : `float`, optional
            History length in number of years used for calibration. A
            fractional number will be rounded to an integer number of months.
            The default is `3.25`.
        pclose : Boolena, optiona \n
            `True` : assumes `mktdata` contains closing prices only, 
            with columns the asset symbols and indexed by the 
            observation dates, \n
            `False` : assumes `mktdata` is in the usual format
            returned by `azapy.mktData` function.
            
        Returns
        -------
        `None`
        """
        if colname is None:
            if self.colname is None:
                raise ValueError("Unspecified price colum (colname).")
        else:
            self.colname = colname
        
        if freq is None:
            if self.freq is None:
                raise ValueError("freq not set - it must be 'Q' or 'M'.")
        elif freq in ['Q', 'M']:
            self.freq = freq
        else:
            raise ValueError(f"wrrong value for freq: {freq}" 
                              " - it must be 'Q' or 'M'.")
            
        if hlength is None:
            if self.hlength is None:
                raise ValueError("hlength must be set to a positive "
                                 "value eq 1")
        else:
            self.hlength = hlength

        if mktdata is None:
            # nothing else to do
            return
        
        periods = 63 if self.freq == 'Q' else 21
        sdate = mktdata.index[int(-np.round(self.hlength * 12) * 21)]

        if pclose == False:
            if self.colname not in mktdata.columns:
                raise ValueError(f"colname: {colname} not in mktdata.columns")
            rrate = mktdata.pivot(columns='symbol', 
                                  values=self.colname)[sdate:]
        else:
            rrate = mktdata[sdate:]
            
        self.last_data = rrate.iloc[-1]
        rrate = rrate.pct_change(periods=periods).dropna()
        self.set_rrate(rrate)


    def getPositions(self, nshares=None, cash=0., ww=None, verbose=True):
        """
        Computes the rebalanced number of shares.

        Parameters
        ----------
        nshares : `panda.Series`, optional
            Initial number of shares per portfolio component.
            A missing component
            entry will be considered 0. A `None` value assumes that all
            components entries are 0. The name of the components must be
            present in the mrkdata. The default is `None`.
        cash : `float`, optional
            Additional cash to be added to the capital. A
            negative entry assumes a reduction in the total capital
            available for rebalance. The total capital cannot be < 0.
            The default is 0. 
        ww : `panda.Series`, optional
            External overwrite portfolio weights. 
            If it not set to `None` these
            weights will overwrite the calibration results.
            The default is `None`. 
        verbose : Boolean, optional
            Is it set to `True` the function prints the closing prices date.
            The default is `True`.
        
        Returns
        -------
        `pandas.DataFrame` : the rolling information.

        Columns:

            - `'old_nsh'` :
                initial number of shares per portfolio component and
                the additional cash. These are input values.
            - `'new_nsh'` :
                the new number of shares per component plus the residual
                cash (due to the rounding to an integer number of shares).
                A negative entry means that the investor needs to add more
                cash to cover for the roundup shortfall.
                It has a small value.
            - `'diff_nsh'` :
                number of shares (buy/sale) needed to rebalance the 
                portfolio.
            - `'weights'` :
                portfolio weights used for rebalancing. The cash entry is
                the new portfolio value (invested capital).
            - `'prices'` :
                the share prices used for rebalances evaluations.

            Note: Since the prices are closing prices, the rebalance can be
            computed after the market close and before the 
            trading execution (next day). 
            Additional cash slippage may occur due
            to share price differential between the previous day closing and
            execution time.
        """
        ns = pd.Series(0, index=self.rrate.columns)
        if nshares is not None:
            ns = ns.add(nshares, fill_value=0)

        if len(self.rrate.columns) != len(ns):
            raise ValueError("nshares must by a subset "
                             f"of {self.rrate.columns}")

        if self.last_data is None:
            raise ValueError("The obj was not set with mktdata")
            
        if verbose:
            print(f"Use closing prices as of {self.last_data.name}")

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


    def _set_getWeights(self, mktdata=None, **params):
        self._reset_output()
        self.verbose = params['verbose'] if 'verbose' in params.keys() else False

        if 'pclose' in params.keys():
            self.set_mktdata(mktdata, pclose=params['pclose'])
        else:
            self.set_rrate(mktdata)


    def _reset_output(self):
        self.status = None
        self.ww = None
        self.time_level1 = None
        self.time_level2 = None
        self.time_level3 = None
        
        
    def getWeights(self):
        pass