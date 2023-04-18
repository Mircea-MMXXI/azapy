import pandas as pd
import numpy as np

from .NullSelector import NullSelector
# from azapy.MkT.MkTcalendar import NYSEgen


class DualMomentumSelector(NullSelector):
    
    def __init__(self, pname='DualMomentum', ftype='13612w', fw=[4, 3, 2, 1], 
                 nw=3, threshold=6, col_price='adjusted'):#, calendar=None):
        super().__init__()
        # need test the input integrity
        self.pname = pname
        self.filter_type = ftype
        self.fw = fw
        self.threshold = threshold
        self.nw = nw
        self.col_price = col_price
        # self.calendar = calendar
        
        # if self.calendar is None:
        #     self.calendar = NYSEgen()
        self.rank = None
        self.mkt = None
        
        
    def getSelection(self, mktdata, **params):
        verbose = params['verbose'] if 'verbose' in params.keys() else False
        
        self.mkt = mktdata.pivot(columns='symbol', values=self.col_price).dropna()
        self.rank = self._filter_rank()
        prank = sum(self.rank > 0)
        capital = min(prank / self.threshold, 1) * min(prank / self.nw, 1)
        pmkt = mktdata.loc[mktdata['symbol'].isin(self.rank.index[:min(prank, self.nw)])]
        
        if verbose: 
            print(f"Selctor {self.pname} :\n\t capital {capital}\n"
                  f"\t selction {pmkt['symbol'].unique()}")
            
        return capital, pmkt
        
    
    def _filter_rank(self):
        if self.filter_type == '13612w':
            filter_func = self._filter_13612w
            
        frank = self.mkt.apply(filter_func, axis=0)
        frank.name = 'rank'
        
        return frank.sort_values(ascending=False)


    def _filter_13612w(self, mktd):
        #edate = mktd.index[-1]
        nmonths = [1, 3, 6, 12]
        ww = self.fw
       
        rrfilter = 0
        for ii, nm in enumerate(nmonths):
            # sdate = edate - pd.offsets.DateOffset(months=nm)
            # sdate = np.busday_offset(sdate.date(), 0, roll='backward',  
            #                          busdaycal=self.calendar)
            sdate = mktd.index[-int(nm * 21)]
            if sdate < mktd.index[0]:
                raise ValueError(f"Not enough data for calibration "
                                 f"{sdate} < first market data record "
                                 f"{mktd.index[0]}")
            rrfilter += ((mktd.iloc[-int(nm * 21)] / mktd.iloc[-1]) ** (12 / nm) - 1) * ww[ii]
        
        return rrfilter / np.sum(ww)
