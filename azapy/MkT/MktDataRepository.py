import pandas as pd 
from copy import deepcopy

from azapy.MkT.MkTcalendar import NYSEgen
from azapy.MkT.MkTreader import MkTreader


class symbolData:
    def __init__(self, name, dates, values,
                 request_time, source, sdate, edate, nr_obs, 
                 error_status, error_log):
        self.name = name
        self.dates = dates
        self.values = values
        
        self.request_time = request_time
        self.source = source
        self.sdate = sdate
        self.edate = edate
        self.error_status = error_status
        self.nr_obs = nr_obs
        
        self.error_log = error_log


class MkTDataRepository:
    def __init__(self):
        self.mktdata = {}
        
        self.source = None
        self.calendar = NYSEgen()
        self._bday = pd.tseries.offsets.CustomBusinessDay(calendar=self.calendar)
        self.file_dir = "MkTData"
        self.file_format = "csv"
        self.api_key = None
        self.param = None
        self.verbose = True
        self.force = False
        self.save = True
        
        
    def get(self, symb, sdate="2012-01-01", edate='today', out_format='frame'):
        rout = {}
        if isinstance(symb, str):
            symb = [symb]
        elif not isinstance(symb, list):
            print(f"Should not be here: Error symbol type: {type(symb)} "
                          + "must be str or a list of str")
        symb_adj = [sy.upper() for sy in symb]
        sdate_adj = pd.to_datetime(sdate).normalize().tz_localize(None)
        sdate_adj = self._bday.rollforward(sdate_adj)
        edate_adj = pd.to_datetime(edate).normalize().tz_localize(None)
        edate_adj = self._bday.rollback(edate_adj)
        
        self._update(symb_adj, sdate_adj, edate_adj)
        
        for sy in symb_adj:
            val = self.mktdata[sy]
            istart = val.dates.get_loc(sdate_adj)
            iend = val.dates.get_loc(edate_adj) + 1
            rout[sy] = deepcopy(val)
            rout[sy].dates = rout[sy].dates[istart:iend]
            rout[sy].values = rout[sy].values.iloc[istart:iend, :]
            
        if out_format in ['frame', 'dict']:
            dict_out = {}
            for sy, val in rout.items():
                dd = val.values
                dd.index = val.dates
                dd.insert(0, 'symbol', sy)
                dict_out[sy] = dd
            if out_format == 'dict':
                return dict_out
            else:
                return pd.concat(dict_out.values())
            
        return rout
            
        
    def _update(self, symb, sdate, edate):
        symb_in = []
        for sy in symb:
            if sy in self.mktdata.keys():
                if ((sdate >= self.mktdata[sy].dates[0]) &
                    (edate <= self.mktdata[sy].dates[-1])):
                    # mktdata for this symbol is in repository
                    symb_in.append(sy)
                
        symb_out = [sy for sy in symb if sy not in symb_in]
        print(f"symb_out {symb_out}")
        mktrd = MkTreader()
        mkt = mktrd.get(symb_out, sdate=sdate, edate=edate, 
                        calendar=self.calendar, output_format='dict',
                        force=self.force, save=self.save,
                        file_format=self.file_format, api_key=self.api_key,
                        param=self.param, verbose=self.verbose)
        mktrep = mktrd.get_request_status()
        error_log = mktrd.get_error_log()
        for kk, val in mkt.items():
           self.mktdata[kk] = symbolData(
                kk, 
                val.index, 
                val.reset_index(drop='True').drop(columns='symbol'),
                request_time = pd.Timestamp.today(),
                source = mktrep.loc['source',kk],
                sdate = mktrep.loc['sdate', kk],
                edate = mktrep.loc['edate', kk],
                nr_obs = mktrep.loc['nrow', kk],
                error_status = mktrep.loc['error', kk],
                error_log = {} if kk not in error_log.keys() else error_log[kk])
                           
       
                            
                            
#==============================================================================
if __name__ == "__main__":
    sdate = "2012-01-01"
    edate = '2024-01-03'
    
    mktr = MkTDataRepository()
    
    symb = ['GLD', 'TLT', 'XLV']
    d1 = mktr.get(symb, sdate, edate, out_format='dict')
                       
    symb2= ['XLV', 'VGT', 'IHI'] 
    d2 = mktr.get(symb2, sdate, edate, out_format='frame')
    
    # import xlwings as xw

    # xw.view(pd.concat(d1.values()))    
    # print(d1)
    
    
