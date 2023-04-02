
class NullSelector:
    def __init__(self):
        self._ptype = 'Selector'
        self.mktdata = None
        self.colname = None
        

    def _set_mktdata(self, mktdata, colname):
        if mktdata is None:
            raise ValueError("No mktdata!!")
        
        if (colname is not None) and (colname not in mktdata.columns):
           raise ValueError(f"colname {colname} not in mktdata.columns")
       
        self.mktdata = mktdata
        self.colname = colname
    
    
    def getSelection(self, mktdata, colname='adjusted'):
        self._set_mktdata(mktdata, colname)
        return 1, self.mktdata