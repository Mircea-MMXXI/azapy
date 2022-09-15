"""
Contains:

    - NYSEgen : NYSE business calendar
"""

import numpy as np
import pandas_market_calendars as mcal

def NYSEgen(sdate='1980-01-01', edate='2050-12-31'):
    """
    Generate `numpy.busdaycalendar` for NYSE between `sdate` and 
    `edate`.

    Parameters
    ----------
    `sdate` : str date, optional
        Calendar start date. The default is `'1980-01-01'`.
    `edate` : str date, optional
        Calendar end date. The default is `'2050-12-31'`.
    Returns
    -------
    `numpy.busdaycalendar`
        NYSE business calendar.
    """
    sdate_ = np.datetime64(sdate)
    edate_ = np.datetime64(edate)
    # get the NYSE holiday list from standard pandas_market_calendars
    ldates = list(mcal.get_calendar('NYSE').holidays().holidays)
    # append here additional holidays 
    # already added to pandas_market_calendars new version
    # ldates.append(np.datetime64('2012-10-29'))
    # ldates.append(np.datetime64('2012-10-30'))
    # ldates.append(np.datetime64('2001-09-11'))
    ldates.sort()
    # move to np.array and return the business calendar
    hdates = np.array(ldates)
    hdates = hdates[(hdates >= sdate_) & ( hdates <= edate_)]
    return np.busdaycalendar(holidays=hdates)