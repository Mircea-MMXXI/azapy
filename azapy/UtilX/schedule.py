""""""  # <- add this
"""
Contains:
    - schedule_simple : creates a simple schedule
    - schedule_roll : creates a rolling history schedule
    - schedule_offset : create a simple schedule with an offset start
"""
import pandas as pd
import numpy as np
import pandas.tseries.offsets as pt

from azapy.MkT.MkTcalendar import NYSEgen


def schedule_simple(sdate='2010-01-01',
                    edate='today',
                    freq='Q',
                    noffset=-3,
                    fixoffset=-1,
                    calendar=None):
    """
    Creates a simple schedule: 'Droll', 'Dfix'.

    Parameters
    ----------
    sdate : `str`, optional
        Start date (reference) of the schedule.
        The default is `'2010-01-01'`.
    edate : `str`, optional
        End date (reference) of the schedule. 
        The default is `'today'`.
    freq : `str`, optional
        Rolling period. It can take 2 values: `'Q'` for quarterly and `'M'` for
        monthly rolling periods. The default is `'Q'`.
    noffset : `int`, optional
        Offset in number of business days for `Droll` relative to the end of 
        calendar period (quarter or month). The default is `-3`.
    fixoffset : `int`, optional
        Offset in number of business days for `Dfix` relative to `Droll`. It 
        must be <=0. The default is `-1`.
    calendar : `numpy.busdaycalendar`, optional
        Business days calendar. If is it `None` then the calendar will be set
        to NYSE business calendar.
        The default is `None`.

    Returns
    -------
    `pandas.DataFrame` : Table containing 2 datetime columns,
        - 'Droll' the rolling date,
        - 'Dfix' the fixing date.
    """
    sdate = pd.to_datetime(sdate)
    if freq == 'Q': edate = pd.to_datetime(edate) + pt.QuarterEnd(1)
    elif freq == 'M': edate = pd.to_datetime(edate) + pt.MonthEnd(1)
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


def schedule_roll(sdate='2010-01-01',
                  edate='today',
                  freq='Q',
                  noffset=-3,
                  fixoffset=-1,
                  calendar=None,
                  hlength=1.25):
    """
    Creates a schedule with rolling history: 'Droll', 'Dfix' and 'Dhist'.

    Parameters
    ----------
    sdate : `str`, optional
        Start date (reference) of the schedule.
        The default is `'2010-01-01'`.
    edate : `str`, optional
        End date (reference) of the schedule. 
        The default is `'today'`.
    freq : `str`, optional
        Rolling period. It can take 2 values: `'Q'` for quarterly and `'M'` for
        monthly rolling periods. The default is `'Q'`.
    noffset : `int`, optional
        Offset in number of business days for Droll relative to the end of 
        calendar period (quarter or month). The default is `-3`.
    fixoffset : `int`, optional
        Offset in number of business days for `Dfix` relative to `Droll`. It 
        can be zero or negative. The default is `-1`.
    calendar : `numpy.busdaycalendar`, optional
        Business days calendar. If is it `None` then the calendar will be set
        to NYSE business calendar.
        The default is `None`.
    hlength : `float`, optional
        Offset in number of years for 'Dhist' relative to 'Dfix'. A fractional 
        value will be rounded to an integer number of months via
        `round(hlength * 12, 0)`. `hlength` must be non-negative.
        The default is `1.25` years.

    Returns
    -------
    `pandas.DataFrame` : Table containing 3 datetime columns, 
        - 'Droll' the rolling date,
        - 'Dfix' the fixing date,
        - 'Dhist' start date for calibration period.
    """
    sdate = pd.to_datetime(sdate)
    edate = pd.to_datetime(edate)
    if calendar is None:
        calendar = NYSEgen()
    sch = schedule_simple(sdate, edate, freq, noffset, fixoffset, calendar)
    sch['Dhist'] = sch['Dfix'] - pd.offsets \
        .DateOffset(months=round(hlength * 12, 0))
    sch['Dhist'] = np.busday_offset(sch['Dhist'].to_numpy(dtype='<M8[D]'), 
                                    0, roll='backward', busdaycal=calendar)
    for k in range(len(sch)):
        if sch['Dhist'][k] >= sdate: 
            schedule = sch[k:].reset_index(drop=True)
            return schedule
    
    raise ValueError("Cannot make a schedule!!")
    
    
def schedule_offset(sdate='2010-01-01',
                    edate='today',
                    freq='Q',
                    noffset=-3,
                    fixoffset=-1,
                    calendar=None,
                    hlength=1.25):
    """
    Creates a simple schedule ('Droll', 'Dfix') with an offset start.

    Parameters
    ----------
    sdate : `str`, optional
        Start date (reference) to which the offset is added.
        The default is `'2010-01-01'`.
    edate : `str`, optional
        End date (reference) of the schedule. 
        The default is `'today'`.
    freq : `str`, optional
        Rolling period. It can take 2 values: `'Q'` for quarterly and `'M'` for
        monthly rolling periods. The default is `'Q'`.
    noffset : `int`, optional
        Offset in number of business days for Droll relative to the end of 
        calendar period (quarter or month). The default is `-3`.
    fixoffset : `int`, optional
        Offset in number of business days for `Dfix` relative to `Droll`. It 
        can be zero or negative. The default is `-1`.
    calendar : `numpy.busdaycalendar`, optional
        Business days calendar. If is it `None` then the calendar will be set
        to NYSE business calendar.
        The default is `None`.
    hlength : `float`, optional
        Offset, in number of years, for first 'Droll' relative to 'sdate'.  
        A fractional value will be rounded to an integer number of months via
        `round(hlength * 12, 0)`. `hlength` must be non-negative.
        The default is `1.25` years.

    Returns
    -------
    `pandas.DataFrame` : Table containing 2 datetime columns
        - 'Droll' the rolling date,
        - 'Dfix' the fixing date.
    """
    sdate = np.datetime64(pd.to_datetime(sdate).date())
    edate = np.datetime64(pd.to_datetime(edate).date())
    if calendar is None:
        calendar = NYSEgen()
    sdate = sdate + pd.offsets.DateOffset(months=round(hlength * 12, 0))
    sdate = np.busday_offset(np.array(sdate, dtype='<M8[D]'), 
                             0, roll='forward', busdaycal=calendar)
    if sdate >= edate:
        raise ValueError(f"offset sdate {sdate} > edate {edate}!!!")
        
    return schedule_simple(sdate, edate, freq, noffset, fixoffset, calendar)
