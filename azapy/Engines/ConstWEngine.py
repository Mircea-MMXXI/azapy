import pandas as pd

class ConstWEngine():
    """
    Constant weighted portfolio.

    **Attributes**
        * status : `int` - computation status (`0` - success, any other 
          value indicates an error)
        * ww : `pandas.Series` - portfolio weights
        * name : `str` - portfolio name
    """
    def __init__(self, ww=None, name=None):
        """
        Constant Weighted Portfolio.

        Parameters
        ----------
        ww : `pandas.Series`
            Portfolio weights (per symbol). If it is `None` the Equal 
            Weighted Portfolio (EWP) will be considered.
            The default is `None`.
        name : `str`, optional
             Portfolio name. The default value is `None`.

        Returns
        -------
        The object.
        """
        self._ptype_ = 'Optimizer'
        self.ww = ww / ww.sum(numeric_only=True)
        self.status = 0
        self.name = name
        
    def getWeights(self, mktdata, **params):
        """
        Returns the portfolio weights adjusted to `mktdata`.\n
        Note: the symbols from `mktdata` take precedents. 
        If `ww` was set to `None` (i.e., EWP) then the weights will be
        set to 1/n, where n is the number symbols in the `mktdata`.
        Otherwise, `ww` will be reduced to the symbols included in `mktdata`
        and renormalized to unit.

        Parameters
        ----------
        mktdata : `pandas.DataFrame`, optional
            The portfolio components historical, prices or rates of return, 
            see`'pclose'` definition below. The default is `None`. \n
            Note: in this call only the values of the symbols are relevant.
        **params: other optional parameters
            `pclose` : Boolean, optional
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
        `pandas.Series` : with portfolio weights.
        """
        if mktdata is None:
            return self.ww
        
        if 'pclose' in params.keys():
            if params['pclose'] == False:
                symb = mktdata['symbol'].unique()
            else:
                symb = mktdata.index
                
        if self.ww is None:
            ww = 1 / len(symb)
            return pd.Series(ww, index=symb)
        
        ww = self.ww[symb]
        return ww / ww.sum(numeric_only=True)
       
        
    def set_mktdata(self, *dummy1, **dummy2):
        """
        Dummy function provided for compatibility with other weights 
        generating (optimizer) classes.

        Returns
        -------
        0
        """
        return self.status
