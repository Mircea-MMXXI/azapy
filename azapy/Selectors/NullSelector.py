
class NullSelector:
    def __init__(self, pname='NullSelector'):
        self._ptype_ = 'Selector'
        self.pname = pname

    
    def getSelection(self, mktdata, **params):
        verbose = params['verbose'] if 'verboss' in params.keys() else False
        if verbose: 
            print(f"Selctor {self.pname} :\n\t capital {1}\n"
                  "\t selction {mktdata['symbols'].unique()}")
        return 1, mktdata