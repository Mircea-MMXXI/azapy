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
    """
    Creates a simple schedule 'Droll', 'Dfix'

    Parameters
    ----------
    sdate : datetime, optional
        Start date (reference) of the schedule.
        The default is pd.to_datetime("2010-01-01").
    edate : datetime, optional
        End date (reference) of the schedule. 
        The default is pd.to_datetime("today").
    freq : string, optional
        Rolling period. It can take 2 values: 'Q' for quarterly and 'M' for
        monthly rolling periods. The default is 'Q'.
    noffset : int, optional
        Offset in number of business days for Droll relative to the end of 
        calendar period (quarter or month). The default is -3.
    fixoffset : int, optional
        Offset in number of business days for Dfix relative to Droll. It 
        can be zero or negative. The default is -1.
    calendar : np.busdaycalendar, optional
        Business days calendar. If is it None then the calendar will be set
        to NYSE business calendar via a call to azapy.NYSEgen(). 
        The default is None.

    Returns
    -------
    pd.DataFrame
        Table containing 2 datetime columns: 'Droll' the rolling date and
        'Dfix' the fixing date.
    """
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
    """
    Creates a schedule with rolling history: 'Droll', 'Dfix' and 'Dhist'

    Parameters
    ----------
    sdate : datetime, optional
        Start date (reference) of the schedule.
        The default is pd.to_datetime("2010-01-01").
    edate : datetime, optional
        End date (reference) of the schedule. 
        The default is pd.to_datetime("today").
    freq : string, optional
        Rolling period. It can take 2 values: 'Q' for quarterly and 'M' for
        monthly rolling periods. The default is 'Q'.
    noffset : int, optional
        Offset in number of business days for Droll relative to the end of 
        calendar period (quarter or month). The default is -3.
    fixoffset : int, optional
        Offset in number of business days for Dfix relative to Droll. It 
        can be zero or negative. The default is -1.
    calendar : np.busdaycalendar, optional
        Business days calendar. If is it None then the calendar will be set
        to NYSE business calendar via a call to azapy.NYSEgen(). 
        The default is None.
    hlength : float, optional
        Offset in number of years for 'Dhist' relative to 'Dfix'. A fractional 
        value will be rounded to an integer number of months via
        round(hlength * 12, 0). hlength must be non negative.
        The default is 1.25.

    Raises
    ------
    ValueError
        The schedule can not be build with the present input parameters. 

    Returns
    -------
    pd.DataFrame
        Table containing 2 datetime columns: 'Droll' the rolling date,
        'Dfix' the fixing date and 'Dhist' start day for a calibration period.
    """
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