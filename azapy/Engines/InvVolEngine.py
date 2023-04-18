import time

from ._RiskEngine import _RiskEngine

class InvVolEngine(_RiskEngine):
    """
    Inverse volatility portfolio.
    
    Methods:
        * getWeights
        * getPositions
        * set_rrate
        * set_mktdata   
    Attributes:
        * status
        * ww
        * name
    """
    def getWeights(self, mktdata=None, **params): 
        """
        Computes the optimal portfolio weights.

        Parameters
        ----------
        `mktdata` : `pandas.DataFrame`, optional;
            The portfolio components historical, prices or rates of return, see
            `'pclose'` definition below.
            If it is not `None`, it will overwrite the set of historical rates
            of return computed in the constructor from `'mktdata'`. 
            The default is `None`. 
        `params`: other optional paramters;
            Most common: \n
            `verbose` : Boolean, optional;
                If it set to `True`, then it will print a messages when 
                the optimal portfolio degenerates to a single asset.
                The default is `False`.
            `pclose` : Boolean, optional;
                If it is absent then the `mktdata` is considered to contain 
                rates of return, with columns the asset symbols and indexed 
                by the observation dates, \n
                `True` : assumes `mktdata` contains closing prices only, 
                with columns the asset symbols and indexed by the 
                observation dates, \n
                `False` : assumes `mktdata` is in the usual format
                returned by `azapy.mktData` function.
  
        Returns
        -------
        `pandas.Series`
            Portfolio weights.
        """
        toc = time.perf_counter()
        self._set_getWeights(mktdata, **params)
        self._calc_ww()
        self.status = 0
        self.time_level1 = time.perf_counter() - toc
        
        return self.ww
    
    
    def _calc_ww(self):
        self.ww = 1. / self.rrate.std(numeric_only=True)
        self.ww /= self.ww.sum(numeric_only=True)