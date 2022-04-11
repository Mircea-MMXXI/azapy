"""
Contains 

- NYSEgen : NYSE business calendar
"""

import numpy as np
import pandas_market_calendars as mcal

def NYSEgen(sdate = np.datetime64('1980-01-01'),
            edate = np.datetime64('2050-12-31')):
    """
    Generate the numpy business calendar from NYSE between ``sdate`` and 
    ``edate``.

    Parameters
    ----------
    sdate : np.datetime64, optional
        Calendar start date. The default is np.datetime64('1980-01-01').
    edate : np.datetime64, optional
        Calendar end date. The default is np.datetime64('2050-12-31').
    Returns
    -------
    numpy.busdaycalendar
        NYSE business calendar.
        It includes as holidays 2001-09-11, 2012-10-29 and 2012-10-30 when
        NYSE was closed.
    """
    # get the NYSE holiday list from standard pandas_market_calendars
    # it doesn't include the days when the exchange was closed
    ldates = list(mcal.get_calendar('NYSE').holidays().holidays)
    # append the exception days when the exchange was closed
    ldates.append(np.datetime64('2012-10-29'))
    ldates.append(np.datetime64('2012-10-30'))
    ldates.append(np.datetime64('2001-09-11'))
    ldates.sort()
    # move to np.array and return the business calendar
    hdates = np.array(ldates)
    hdates = hdates[(hdates >= sdate) & ( hdates <= edate)]
    return np.busdaycalendar(holidays=hdates)