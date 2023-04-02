import pandas as pd
import numpy as np

from .NullSelector import NullSelector
from azapy.MkT.MkTcalendar import NYSEgen


class DualMomentumSelector(NullSelector):
    
    def __init__(self, ftype='13612w', fw=[4, 3, 2, 1], nw=3, threshold=6, 
                 calendar=None):
        super().__init__()
        # need test the input integrity
        self.filter_type = ftype
        self.fw = fw
        self.threshold = threshold
        self.nw = nw
        
        if calendar is None:
            self.calendar = NYSEgen()
        else:
            self.calendar = calendar
        
        self.rank = None
        self.mkt = None
        

 
    def _filter_13612w(self, mktd):
        edate = mktd.index[-1]
        #ppe = mktd[edate]
        nmonths = [1, 3, 6, 12]
        ww = self.fw
       
        rrfilter = 0
        for ii, nm in enumerate(nmonths):
            sdate = edate - pd.offsets.DateOffset(months=nm)
            sdate = np.busday_offset(sdate.date(), 0, roll='backward',  
                                     busdaycal=self.calendar)
            #pps = mktd[sdate]
            rrfilter += ((mktd[edate] / mktd[sdate]) ** (12 / nm) - 1) * ww[ii]
        
        return rrfilter / np.sum(ww)
    
    def _filter_rank(self):
        if self.filter_type == '13612w':
            filter_func = self._filter_13612w
            
        frank = self.mkt.apply(filter_func, axis=0)
        frank.name = 'rank'
        
        return frank.sort_values(ascending=False)
    
    def getSelection(self, mktdata=None, colname='adjusted'):
        self._set_mktdata(mktdata, colname)
        if self.colname is None:
            self.mkt = self.mktdata.copy()
        else:
            self.mkt = self.mktdata.pivot(columns='symbol', values=self.colname).dropna()
        self.rank = self._filter_rank()
        prank = sum(self.rank > 0)
        capital = min(prank / self.threshold, 1) * min(prank / self.nw, 1)
        if self.colname is None:
            pmkt = self.mktdata[self.rank.index[:min(prank, self.nw)]]
        else: 
            pmkt = self.mktdata.loc[self.mktdata['symbol'].isin(self.rank.index[:min(prank, self.nw)])]
    
        return capital, pmkt

    
        