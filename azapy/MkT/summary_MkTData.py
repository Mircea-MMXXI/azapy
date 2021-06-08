# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 17:35:22 2021

@author: mircea
"""

import pandas as pd
import numpy as np
from collections import defaultdict 

from .readMkTData import NYSEgen

def summary_MkTData(rprice, calendar=None, sdate=None, edate=None):
    """
    Summary of MkT Dats time-series lentgh and quality

    Parameters
    ----------
    rprice : dict or pd.DataFrame 
        MkT Data in the format returned by azapy.readMkT function.
    calendar : np.busdaycalendar, optional
        Business days calendar. If is set to None it will 
        default to NYSE business calendar.
    sdate : pd.Timestamp, optional
        Time-series start date. If is None then sdate will be set to the 
        earliest date in rprice.
        Default value is None.
    edate : pd.Timestamp, optional
        Time-earies end date. If it None then edate will be set to the most 
        recent date in rprice.
        Default value is None.

    Returns
    -------
    TYPE pd.DataFrame
        Tabel with columns:
            symbol - time-series symbol
            begin - start date
            end - end date
            length - number of records
            na - total number of NA
            na_b - number of missing records at the begining
            na_e = number of missing records at the end
            cont - total number of missing records
            
    
    Comment: the main aplication is to asses the missing data in the 
    time-series extracted with azapy.readMkT function.
    """
    
    if isinstance(rprice, dict):
        gite = rprice.items()
    else:
        gite = rprice.groupby('symbol')
    
    sds = []
    eds = []
    for _ , v in gite:
        sds.append(v.index[0])
        eds.append(v.index[-1])
    msds = min(sds)
    meds = max(eds)
    if (sdate is None) or (sdate > msds): sdate = msds
    if (edate is None) or (edate < meds): edate = meds
    
    if calendar is None: calendar = NYSEgen()

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
        res['na'].append(v.isnull().sum().sum())
        res['na_b'].append(0 if sd == hd[0] else hd.get_loc(sd) )
        res['na_e'].append(0 if ed == hd[-1] else len(hd) - hd.get_loc(ed) - 1)
        res['cont'].append(len(hd) - len(v.index))

    return pd.DataFrame(res)
