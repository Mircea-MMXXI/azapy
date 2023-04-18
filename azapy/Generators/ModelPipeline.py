import pandas as pd
import numpy as np

from azapy.Selectors.NullSelector import NullSelector

class ModelPipeline:
    """
    Construct a portfolio model from a sequnce os elementary 
    models. The last element of the squence must be an optomizer model 
    while the rest could be any nymber of selector models. 
    """
    def __init__(self, sequence=[NullSelector(), "EWP"]):
        """
        Constructor

        Parameters
        ----------
        sequence : list, optional;
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
            raise ValueError("sequence cannot be emply!")
        self.sequence = sequence
        self.assets_opt = None
            
        
    def getWeights(self, mktdata, **params):
        """
        Computes the portfolio weights.

        Parameters
        ----------
        mktdata : `pandas.DataFrame`;
            Historical daily market data as returned by `azapy.readMkT` 
            function.
        **params : 
            Additional parameters that may be requierd by the elementatry 
            models. An example is `verbose=True`.

        Returns
        -------
        ww : `pandas.Series`:
            Portfolio weights.
        """
        if mktdata is None:
            raise ValueError("no mktdata!")
        mktd = mktdata
        
        verbose = params['verbose'] if 'verbose' in params.keys() else False
 
        anames = np.append(mktd['symbol'].unique(),'_CASH_')
        ww = pd.Series(np.zeros(len(anames)), index=anames)
        
        # loop through selectrs
        capital = 1
        for k in range(self.lseq - 1):
            if self.sequence[k]._ptype_ != 'Selector':
                raise ValueError("unknown obj as a `Selector` !!!")
                
            cpt, mktd = self.sequence[k].getSelection(mktd, **params)
            capital *= cpt
            if capital == 0:
                # all in cash
                ww['_CASH_'] = 1
                return ww

        # selected assets
        anames = np.append(mktd['symbol'].unique(),'_CASH_')
        nassets = len(anames)
            
        self.assets_opt = anames
        if verbose:
            print(f"assets for optimization: {self.assets_opt[:-1]}"
                  f"\nwith calpital in cash {1-capital}")
        
        # allocation strategy
        
        # strategies by name (str)
        if type(self.sequence[-1]) == str:
            # EWP (equal weighted portfolio)
            if self.sequence[-1] == 'EWP':
                ww[anames[:-1]] = capital / (nassets - 1)
                ww['_CASH_'] = 1 - capital
                return ww
            else:
                raise ValueError(f"unknown optimizer {self.sequence[-1]}")
        # Analyzer (optimizer)
        elif self.sequence[-1]._ptype_ == 'Optimizer':
            params['pclose'] = False
            ws = self.sequence[-1].getWeights(mktdata=mktd, **params)
            ww[ws.index] = ws * capital
            ww['_CASH_'] = 1 - capital
            return ww
        else:
            raise ValueError("unknown optimizer")
            