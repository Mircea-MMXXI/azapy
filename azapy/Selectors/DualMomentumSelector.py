import numpy as np

from .NullSelector import NullSelector


class DualMomentumSelector(NullSelector):
    """
    Dual Momentum Selector.
    
    Given a filter it selects the best candidates among the once with 
    highest moment. 
    
    **Attributes**
        * pname : `str` - portfolio name
        * mkt : `pandas.DataFrame` - selection's market data
        * rank : `pandas.Series` - filter rank of all symbols 
        * symb : `list` - selected symbols
        * symb_omitted : `list` - unselected symbols
        * capital : `float` - capital at risk as a fraction of the total 
          capital
    """
    def __init__(self, pname='DualMomentum', ftype='f13612w', fw=None, 
                 nw=3, threshold=6, col_price='adjusted', **kwargs):
        """
        Constructor

        Parameters
        ----------
        pname : `str`, optional
            Selector name. The default is 'DualMomentum'.
        ftype : `str`, optional
            The filter name (at this point only `'f13612w'` filter is supported). 
            The default is `'f13612w'` are equal weights.
        fw : `list`, optional
            List of filter wights. 
            For `'f13612w'` it must be a list of 4 positive (not all zero)
            numbers. A value of `None` indicates equal weights.
            Note: the weights are normalized internally.
            The default is `None`.
        nw : `int`, optional
            Maximum number of selected symbols. The default is 3.
        threshold : `int`, optional
            Minimum number of symbols with positive momentum for a full 
            capital allocation. The default is 6.
        col_price : `str`, optional
            Name of the price column in the `mktdata` to be used in the 
            momentum evaluations. The default is 'adjusted'.
        **kwargs : `dict`, optional
            Holder for other args.

        Returns
        -------
        The object.
        """
        super().__init__(pname)
        # need test the input integrity
        self.filter_type = ftype
        self.fw = fw
        self.threshold = threshold
        self.nw = nw
        self.col_price = col_price
        self.rank = None
        self._mkt = None
        self.symb = None
        self.symb_omitted = None
        self.capital = None
        
        
    def getSelection(self, mktdata, **params):
        """
        Computes the selection.

        Parameters
        ----------
        mktdata : `pandas.DataFrame`
            MkT data in the format produced by the `azapy` function `readMkT`.
        **params : `dict`, optional
            Other optional parameters:
                **verbose** : `Boolean`, optional
                    When it is set to `True`, the selection symbols are printed.
                    The default is `False`.

        Returns
        -------
        (capital, mkt) : tuple
            capital : `float`
                Fraction of capital allocated to the selection (a positive number 
                <= 1, 1 being full allocation). One minus this value is the
                fraction of reserved capital invested in cash.
            mkt  : `pandas.DataFrame`
                Selection MkT data in the format produced by the `azapy` 
                function `readMkT`.
        """
        verbose = params['verbose'] if 'verbose' in params.keys() else False
        
        self._mkt = mktdata.pivot(columns='symbol', values=self.col_price).dropna()
        self.rank = self._filter_rank()
        prank = sum(self.rank > 0)
        self.symb = self.rank.index[:min(prank, self.nw)].to_list()
        self.capital = min(prank / self.threshold, 1) * min(prank / self.nw, 1)
        self.mkt = mktdata.loc[mktdata['symbol'].isin(self.symb)]
        self.symb_omitted = list(np.setdiff1d(self._mkt.columns, self.symb))
        
        if verbose: 
            print(f"Selector {self.pname} :\n\t capital {self.capital}\n"
                  f"\t selection {self.symb}")
            
        return self.capital, self.mkt
        
    
    def _filter_rank(self):
        if self.filter_type == 'f13612w':
            filter_func = self._f13612w
        else:
            raise ValueError(f"Unknown filter type {self.filter_type}")
            
        frank = self._mkt.apply(filter_func, axis=0)
        frank.name = 'rank'
        
        return frank.sort_values(ascending=False)


    def _f13612w(self, mktd):
        nmonths = [1, 3, 6, 12]
        if self.fw is None:
            ww = [0.25] * 4
        else:    
            ww = self.fw
       
        rrfilter = 0
        for ii, nm in enumerate(nmonths):
            sdate = mktd.index[-int(nm * 21)]
            if sdate < mktd.index[0]:
                raise ValueError(f"Not enough data for calibration "
                                 f"{sdate} < first market data record "
                                 f"{mktd.index[0]}")
            rrfilter += ((mktd.iloc[-int(nm * 21)] 
                          / mktd.iloc[-1]) ** (12 / nm) - 1) * ww[ii]
        
        return rrfilter / np.sum(ww)
