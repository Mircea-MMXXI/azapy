import pandas as pd
import numpy as np
from collections import defaultdict 

from .MkTcalendar import _set_calendar_exchange


def summary_MkTdata(mktdata, calendar=None, sdate=None, edate=None):
    """
    Summary of MkT data time-series length and quality (checks for missing
    records).
    
    Notes
    -----
    Its main application is to assess the missing data in the 
    time-series extracted with `azapy.readMkT` function.

    Parameters
    ----------
    mktdata : `pandas.DataFrame` or a dict of `pandas.DataFrame`
        Market Data in the format returned by `azapy.readMkT` function.
    calendar : `str` or `numpy.busdaycalendar`, optional
            Business calendar. It can be the exchange calendar name as a `str` or 
            a `numpy.busdaycalendar` object.
            If it is `None` then it will be set to NYSE
            business calendar. The default
            value is `None`.
    sdate : date like, optional
        Time-series start date. If it is `None` then `sdate` will be set to the 
        earliest date in mktdata.
        The default is `None`.
    edate : date like, optional
        Time-series end date. If it is `None` then `edate` will be set to 
        the most recent date in mktdata.
        The default is `None`.

    Returns
    -------
    `pandas.DataFrame` : A table with columns:
        - `symbol` : time-series symbol
        - `begin` : start date
        - `end` : end date
        - `length` : number of records
        - `na_total` : total number of `nan`
        - `na_b` : number of missing records at the beginning
        - `na_e` : number of missing records at the end
        - `cont` : total number of missing records
    """
    gite =  mktdata.items() if isinstance(mktdata, dict) else mktdata.groupby('symbol')
    sdate = np.datetime64(sdate) if sdate is not None else min([np.datetime64(v.index[0]) for _, v in gite])
    edate = np.datetime64(edate) if edate is not None else max([np.datetime64(v.index[-1]) for _, v in gite])
    
    calendar = _set_calendar_exchange(calendar)
    hd = pd.DatetimeIndex(pd.bdate_range(sdate, edate, freq='C', holidays=calendar.holidays))
    
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
