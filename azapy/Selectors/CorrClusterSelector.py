import pandas as pd
import numpy as np
from .NullSelector import NullSelector

from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as shc

class CorrClusterSelector(NullSelector):
    """
    Selects symbols with lower inter-correlation.
    
    **Attributes**
        * pname : `str` - portfolio name
        * mkt : `pandas.DataFrame` - selection's market data
        * symb : `list` - selected symbols
        * symb_omitted : `list` - unselected symbols
        * capital : `float` - always set to `1`
    """
    def __init__(self, pname='CorrCluster', corr_threshold=0.95, freq='Q',
                 ftype='f13612w', fw=None, col_price='adjusted', hlength=1):
        """
        Constructor

        Parameters
        ----------
        pname : `str`, optional
            Selector name. The default is 'DualMomentum'.
        corr_thresold : `float`, optional
            Cluster correlation threshold (i.e., a cluster contains only symbols 
            with inter-correlation higher than `corr_threshold`. 
            The default is 0.95.
        freq : `str`, optional
            The horizon of rates subject to correlation estimations. 
            It can be either `'M'` for monthly or `'Q'` for quarterly rates. 
            The defualt is `'Q'`.
        ftype : `str`, optional
            Inner-cluster filter (i.e., criteria to designate the representative
            of a cluster with more than one symbol). At this point only 
            `'f13612w'` is implemented.
            The default is 'f13612w'.
        fw : `list`, optional
            List of filter wights. 
            For `'f13612w'` it must be a list of 4 positive (not all zero)
            numbers. A value of `None` indicates equal weights.
            Note: the weights are normalized internally.
            The default is `None`.
        col_price : `str`, optional
            The name of the pricing column to be considered in computations.
            The default is 'adjusted'.
        hlength : 'float', optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is `3.25` years.
            
        Returns
        -------
        The object.
        """
        super().__init__(pname)
        self.corr_threshold = corr_threshold
        self.freq = freq
        self.filter_name = ftype
        self.fw = fw
        self.col_price = col_price
        self.hlength = hlength
        self.capital = 1 
        self.mkt = None
        self.symb = None
        self.symb_omitted = None
        
        
    def getSelection(self, mktdata, **params):
        """
        Computes the selection.

        Parameters
        ----------
        mktdata : `pandas.DataFrme`
            MkT data in the format produced by the `azapy` function `readMkT`.
        **params : `dict`, optional
            Other optional parameters:
                **verbose** : Boolean, optional
                    When it is set to `True`, the selection symbols are printed.
                    The default is 'False'.
                
                **view** : Boolean, optional
                    If set to `True`, then the dendrogram of hierarchical 
                    classification is printed out. The default is `False`.
                    Note: the tree cutoff is at `1 - corr_threshold` level.

        Returns
        -------
        (capital, mkt) : tuple
            caplital : `float`
                Fraction of capital allocated to the selection. For this 
                selector it is always 1.
            mkt  : `pandas.DataFrame`
               Selection MkT data in the format produced by the `azapy` 
               function `readMkT`.
        """
        mkt = mktdata.pivot(columns='symbol', values=self.col_price)
        symb = mkt.columns
        if self.freq == 'M':
            periods = 21
        elif self.freq == 'Q':
            periods = 63
        else:
            raise ValueError(f"unknown freq {self.freq} -it must be 'M' or 'Q'")
        
        if self.hlength is None:
            raise ValueError("hlength must be set to a positive "
                             "value - e.g. 1")
        idx_sdate = int(-np.round(self.hlength * 12) * 21)
        
        rmkt = mkt.pct_change(periods=periods).dropna()
        scorr = pdist(rmkt.iloc[idx_sdate:].T, metric='correlation')
        Z = shc.linkage(scorr, method='ward', optimal_ordering=True)
        
        view = params['view'] if 'view' in params.keys() else False
        if view:
            _ = shc.dendrogram(Z, labels=symb)
            
        height = 1 - self.corr_threshold
        ctt = shc.cut_tree(Z, height=height)
        sctt = pd.Series(np.transpose(ctt)[0], index=symb)
        
        if self.filter_name == 'f13612w':
            filter_func = self._f13612w
            
        self.symb = []
        for k in range(sctt.max()):
            csy = sctt[sctt == k].index.to_list()
            if len(csy) == 1:
                self.symb.append(csy[0])
            else:
                rank = filter_func(mkt[csy]) 
                self.symb.append(rank.idxmax())

        verbose = params['verbose'] if 'verbose' in params.keys() else False
        if verbose: 
            print(f"Selctor {self.pname} :\n\t capital {1}\n"
                  f"\t selction {self.symb}")
            
        self.mkt = mktdata[mktdata['symbol'].isin(self.symb)]
        self.symb_omitted = list(np.setdiff1d(mkt.columns, self.symb))
                
        return self.capital, self.mkt
         
        
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

