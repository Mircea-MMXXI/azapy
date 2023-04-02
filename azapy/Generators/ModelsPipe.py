import pandas as pd
import numpy as np

from azapy.Selectors.NullSelector import NullSelector

class ModelsPipe:
    def __init__(self, sequence=[NullSelector(), "EWP"]):
        self.lseq =  len(sequence)
        if self.lseq < 1:
            raise ValueError("sequence cannot be emply!")
        self.sequence = sequence
        self.assets_opt = None
            
    def getWeights(self, mktdata, colname='adjusted',
                   # mu=None, d=1, mu0=0., aversion=None,
                   # ww0=None, rrate=None, 
                   verbose=False):
        if mktdata is None:
            raise ValueError("no mktdata!")
        mktd = mktdata
 
        if colname is None:
            anames = mktd.columns.to_list() + ['_CASH_']
        elif colname in mktd.columns:
            anames = np.append(mktd['symbol'].unique(),'_CASH_')
        else:
            raise ValueError("colname is not None or in the mktdata.columns")
        ww = pd.Series(np.zeros(len(anames)), index=anames)
        
        # loopp through selectrs
        capital = 1
        for k in range(self.lseq - 1):
            cpt, mktd = self.sequence[k].getSelection(mktd, colname=colname)
            capital *= cpt
            if capital == 0:
                break
        
        # exit if empty selection
        if capital == 0:
            ww['_CASH_'] = 1
            return ww
        
        # selected assets
        if colname is None:
            anames = mktd.columns.to_list() + ['_CASH_']
        else:
            anames = np.append(mktd['symbol'].unique(),'_CASH_')
        nassets = len(anames)
            
        self.assets_opt = anames
        if verbose:
            print(f"assets for optimization: {self.assets_opt}")
        
        # allocation strategy
        
        # strategies by name (str)
        if type(self.sequence[-1]) == str:
            # EWP (equal weighted portfolio)
            if self.sequence[-1] == 'EWP':
                ww[anames[:-1]] = capital / (nassets - 1)
                ww['_CASH_'] = 1 - capital
                return ww
            
            raise ValueError(f"unknown optimizer {self.sequence[-1]}")
        # Analyzer (optimizer)
        elif self.sequence[-1]._ptype_ == 'Analyzer':
            self.sequence[-1].set_mktdata(mktd, colname)
            # ws = self.sequence[-1].getWeights(mu=mu, mu0=mu0,
            #         aversion=aversion, ww0=ww0, verbose=verbose)
            ws = self.sequence[-1].getWeights(verbose=verbose)
            ww[ws.index] = ws * capital
            ww['_CASH_'] = 1 - capital
            return ww
        else:
            ValueError("unknown optimizer")
            