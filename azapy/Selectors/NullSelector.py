
class NullSelector:
    """
    Null Selector - produces no selection.
    
    Its main purpose is to serve as a base class and to provide consistency 
    across the selector family.
    """
    def __init__(self, pname='NullSelector'):
        """
        Constructor

        Parameters
        ----------
        pname : `str`, optional
            Selector name. The default is 'NullSelector'.

        Returns
        -------
        The object

        """
        self._ptype_ = 'Selector'
        self.pname = pname

    
    def getSelection(self, mktdata, **params):
        """
        Produces the selection.

        Parameters
        ----------
        mktdata : `pandas.DataFrames`
            Market data in the format produced by the `azapy` function 
            `readMKT`.
        **params : `dict`, optional
            Additional optional parameters:
                **verbose** : Boolean, optional
                  When it is set to `True`, the selection symbols are printed.

        Returns
        -------
        (capital, mkt) : tuple
            * capital : `float`, always set to `1`.
            * mkt : `pandas.DataFrame`, always the input `mktdata`
        """
        verbose = params['verbose'] if 'verboss' in params.keys() else False
        if verbose: 
            print(f"Selctor {self.pname} :\n\t capital {1}\n"
                  "\t selction {mktdata['symbols'].unique()}")
        return 1, mktdata