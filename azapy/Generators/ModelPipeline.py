import pandas as pd
import numpy as np

from azapy.Selectors.NullSelector import NullSelector

class ModelPipeline:
    """
    Construct a portfolio weights model from a sequence of elementary models. 
    The last element of the sequence must be an optimizer model while the 
    rest could be any number of selector models. 

    **Attributes**
        - sequence : `list` -  the sequence of elementary models
        - capital : `flat` -  the capital at risk as a fraction of the 
          total capital
        - mktdata : `pandas.DataFrame` - historical market data of selected 
          symbols
        - active_symb : `list` - the list of selected symbols
        - ww : `pandas.Series` - portfolio weights per symbol (all symbols)
        
    **Note**
    All unselected symbols have 0 weight. However, some of the selected symbols
    can have 0 weight after portfolio optimization stage.
    """
    def __init__(self, sequence=[NullSelector(), "EWP"]):
        """
        Constructor

        Parameters
        ----------
        sequence : `list`, optional
            List of elementary models. The last element of the list must be 
            an optimizer while the rest could be any number of selectors. 
            The sequence is executed from right to left.
            The default is `[NullSelector(), "EWP"]`.

        Returns
        -------
        The object.
        """
        self.lseq =  len(sequence)
        if self.lseq < 1:
            raise ValueError("sequence cannot be empty!")
        self.sequence = sequence
        self.assets_opt = None
        self.ww = None
        self.capital = None
        self.mktdata = None
        self.active_symb = None
        self.status = 1
        
        
    def _build_mktdata(self, mktdata):        
        if isinstance(mktdata, dict):
            lvv = []
            for kk, vv in mktdata.items():
                if 'symbol' in vv.columns:
                    lv = vv
                else:
                    lv = vv.copy()
                    lv['symbol'] = kk
                lvv.append(lv)
            return pd.concat(lvv)
        else:
            return mktdata.copy()
            
        
    def getWeights(self, mktdata, **params):
        """
        Computes the portfolio weights.

        Parameters
        ----------
        mktdata : `pandas.DataFrame`
            Historical daily market data as returned by `azapy.readMkT` 
            function.
        **params : optional
            Additional parameters that may be required by the elementary 
            models. An example is `verbose=True`.

        Returns
        -------
        `pandas.Series` : Portfolio weights per symbol.
        """
        if mktdata is None:
            raise ValueError("no mktdata!")
        self.mktdata = self._build_mktdata(mktdata)
        mktd = self.mktdata
        
        verbose = params['verbose'] if 'verbose' in params.keys() else False
 
        anames = np.append(mktd['symbol'].unique(), '_CASH_')
        self.ww = pd.Series(np.zeros(len(anames)), index=anames)
        
        # loop through selectors
        self.capital = 1
        for k in range(self.lseq - 1):
            if pd.isnull(self.sequence[k]):
                continue
            if self.sequence[k]._ptype_ != 'Selector':
                raise ValueError("unknown obj as a `Selector` !!!")
                
            cpt, mktd = self.sequence[k].getSelection(mktd, **params)
            self.capital *= cpt
            if self.capital == 0:
                # all in cash
                self.status = 0
                self.ww['_CASH_'] = 1
                self.active_symb = ['_CASH_']
                return self.ww

        # selected assets
        anames = np.append(mktd['symbol'].unique(),'_CASH_')
        nassets = len(anames)
            
        self.assets_opt = anames
        if verbose:
            print(f"assets for optimization: {self.assets_opt[:-1]}"
                  f"\nwith capital in cash {1-self.capital}")
        
        # allocation strategy
        
        # strategies by name (str)
        if isinstance(self.sequence[-1], str):
            # EWP (equal weighted portfolio)
            if self.sequence[-1] == 'EWP':
                self.ww[anames[:-1]] = self.capital / (nassets - 1)
                self.ww['_CASH_'] = 1 - self.capital
                self.active_symb = self.ww.drop('_CASH_').index.to_list()
                self.status = 0
                return self.ww
            else:
                raise ValueError(f"unknown optimizer {self.sequence[-1]}")
        # Analyzer (optimizer)
        elif self.sequence[-1]._ptype_ == 'Optimizer':
            params['pclose'] = False
            ws = self.sequence[-1].getWeights(mktdata=mktd, **params)
            self.status = self.sequence[-1].status
            if self.status != 0:
                self.ww = None
                self.active_symb = []
            else:
                self.ww[ws.index] = ws * self.capital
                self.ww['_CASH_'] = 1 - self.capital
                self.active_symb = ws.index.to_list()
            return self.ww
        else:
            raise ValueError("unknown optimizer")
            
            
    def getPositions(self, nshares=None, cash=0., ww=None, nsh_round=True, 
                     verbose=True):
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
        nsh_round : `Boolean`, optional
            If it is `True` the invested numbers of shares are round to the 
            nearest integer and the residual cash capital 
            (positive or negative) is carried to the next reinvestment cycle. 
            A value of `False` assumes investments with fractional number 
            of shares (no rounding). The default is `True`.
        verbose : `Boolean`, optional
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
                the share prices used for rebalance evaluations.

            Note: Since the prices are closing prices, the rebalance can be
            computed after the market close and before the 
            trading execution (next day). 
            Additional cash slippage may occur due
            to share price differential between the previous day closing and
            execution time.
        """
        if self.status != 0:
            if verbose:
                print(f"No weights - status {self.status}")
            return None
        
        shares_round = 0 if nsh_round else 16
        price = self.mktdata.pivot(columns='symbol', values='close').iloc[-1]
        price['_CASH_'] = 1            
      
        if ww is None: 
            ww = self.ww
        else:
            ww = ww / ww.sum()
            ww = pd.Series(0, index=price.index).add(ww, fill_value=0)
            
        nsh_in = pd.Series(0, index=price.index)
        if nshares is not None:
            nsh_in = nsh_in.add(nshares, fill_value=0)
            
        cap = nsh_in @ price + cash
        nsh_out = (ww * cap / price)
        nsh_out[nsh_out.index != '_CASH_'] = \
            nsh_out[nsh_out.index != '_CASH_'].round(shares_round)

        cap_adj = cap - nsh_out @ price
 
        nsh_in['_CASH_'] = cash
        nsh_out['_CASH_'] += cap_adj
        
        res = pd.DataFrame({'old_nsh': nsh_in.round(2),
                            'new_nsh': nsh_out.round(2),
                            'diff_nsh': (nsh_out - nsh_in).round(2),
                            'weights' : ww.round(4),
                            'prices': price.round(2)})
        
        res = res.rename(index={'_CASH_': 'cash'})
        res = res[(res['old_nsh'] != 0) | (res['new_nsh'] != 0) | (res.index == 'cash')]
        
        if verbose:
            print(f"Computed as of: {price.name.date()}")
            print(f"previous invested capital: {np.round(cap - cash, 2)}")
            print(f"cash adjustment: {np.round(cash, 2)}")
            print(f"new invested capital: {np.round(cap * self.capital - cap_adj, 2)}")
            print(f"cash in hands: {np.round(cap * (1 - self.capital) + cap_adj, 2)}")
        
        return res