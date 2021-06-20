# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 15:30:50 2021

@author: mircea
"""
import pandas as pd
import numpy as np
import pandas.tseries.offsets as pt

from azapy.MkT.readMkTData import NYSEgen

def simple_schedule(sdate=pd.to_datetime("2010-01-01"),
                    edate=pd.to_datetime("today"),
                    freq='Q',
                    noffset=-3,
                    fixoffset=-1,
                    calendar=None):
    if freq == 'Q': edate = edate + pt.QuarterEnd(1)
    elif freq == 'M': edate = edate + pt.MonthEnd(1)
    else: raise ValueError("Wrong freq, Must be 'Q' or 'M'")
    
    if calendar is None:
        calendar = NYSEgen()
    
    tedx = pd.date_range(start=sdate, end=edate, freq=freq)\
        .to_numpy(dtype='<M8[D]')
    troll = np.busday_offset(tedx, noffset, roll='backward', 
                             busdaycal=calendar)
    tfix = np.busday_offset(troll, fixoffset, roll='backward', 
                            busdaycal=calendar)
    return pd.DataFrame({'Droll': troll, 'Dfix': tfix})

def schedule_roll(sdate=pd.to_datetime("2010-01-01"),
                  edate=pd.to_datetime("today"),
                  freq='Q',
                  noffset=-3,
                  fixoffset=-1,
                  calendar=None,
                  hlength=1.25):
    if calendar is None:
        calendar = NYSEgen()
    sch = simple_schedule(sdate, edate, freq, noffset, fixoffset, calendar)
    sch['Dhist'] = sch['Dfix'] - pd.offsets \
        .DateOffset(months=round(hlength * 12, 0))
    sch['Dhist'] = np.busday_offset(sch['Dhist'].to_numpy(dtype='<M8[D]'), 
                                    0, roll='backward', 
                                    busdaycal=calendar)
    for k in range(len(sch)):
        if sch['Dhist'][k] >= sdate: 
            schedule = sch[k:].reset_index(drop=True)
            return schedule
    
    raise ValueError("Cannot make a schedule!!")