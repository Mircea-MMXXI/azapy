import time

from ._RiskEngine import _RiskEngine

class InvVolEngine(_RiskEngine):
    """
    Inverse volatility portfolio.
    
    **Attributes**
        * status : `int` - computation status (`0` - success, any other 
          value indicates an error)
        * ww : `pandas.Series` - portfolio weights
        * name : `str` - portfolio name
    """
    def getWeights(self, mktdata=None, **params): 
        """
        Computes the optimal portfolio weights.

        Parameters
        ----------
        mktdata : `pandas.DataFrame`, optional
            The portfolio components historical, prices or rates of return, see
            `'pclose'` definition below.
            If it is not `None`, it will overwrite the set of historical rates
            of return computed in the constructor from `'mktdata'`. 
            The default is `None`. 
        **params : other optional parameters
            Most common: \n
            `verbose` : `Boolean`, optional
                If it is set to `True`, then it will print a computation
                messages. The default is `False`.
            `pclose` : `Boolean`, optional
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
        `pandas.Series` : Portfolio weights.
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