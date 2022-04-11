import pandas as pd
import numpy as np
from collections import defaultdict 

#from .readMkTData import NYSEgen
from .MkTcalendar import NYSEgen

def summary_MkTData(mktdata, calendar=None, sdate=None, edate=None):
    """
    Summary of MkT data time-series length and quality (cheks for missing
    records).

    Parameters
    ----------
    mktdata :pd.DataFrame or a dictonary of pd.DataFrame's
        Market Data in the format returned by azapy.readMkT function.
    calendar : np.busdaycalendar, optional
        Business days calendar. If is set to None it will 
        default to NYSE business calendar.
    sdate : pd.Timestamp, optional
        Time-series start date. If it is None then sdate will be set to the 
        earliest date in mktdata.
        The default is None.
    edate : pd.Timestamp, optional
        Time-series end date. If it is None then edate will be set to the most 
        recent date in mktdata.
        The default is None.

    Returns
    -------
    pd.DataFrame: a table with columns:
        - symbol : time-series symbol
        - begin : start date
        - end : end date
        - length : number of records
        - na_total : total number of `nan`
        - na_b : number of missing records at the beginning
        - na_e : number of missing records at the end
        - cont : total number of missing records
            
    Comment: the main application is to asses the missing data in the 
    time-series extracted with azapy.readMkT function.
    """
    if isinstance(mktdata, dict):
        gite = mktdata.items()
    else:
        gite = mktdata.groupby('symbol')
    
    sds = []
    eds = []
    for _ , v in gite:
        sds.append(v.index[0])
        eds.append(v.index[-1])
    msds = min(sds)
    meds = max(eds)
    if (sdate is None) or (sdate > msds): sdate = msds
    if (edate is None) or (edate < meds): edate = meds
    
    if calendar is None: 
        calendar = NYSEgen()

    hd = pd.DatetimeIndex([dd  for dd in pd.date_range(sdate, edate) \
                           if np.is_busday(dd.date(), busdaycal=calendar)])
    
    res = defaultdict(lambda: [])
    for k, v in gite:
        sd = v.index[0]
        ed = v.index[-1]
        res['symbol'].append(k)
        res['begin'].append(sd)
        res['end'].append(ed)
        res['length'].append(len(v))
        res['na_total'].append(v.isnull().sum().sum())
        res['na_b'].append(0 if sd == hd[0] else hd.get_loc(sd) )
        res['na_e'].append(0 if ed == hd[-1] else len(hd) - hd.get_loc(ed) - 1)
        res['cont'].append(len(hd) - len(v.index))

    return pd.DataFrame(res)
