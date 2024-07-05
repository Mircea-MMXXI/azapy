import numpy as np
import pandas_market_calendars as mcal
import warnings


def NYSEgen(sdate='1980-01-01', edate='2050-12-31'):
    """
    Returns the NYSE business calendar between `sdate` and 
    `edate`. 

    `To be deprecated in future versions.
    Instead, please use
    calendarGen(name="NYSE", sdate=sdate, edate=edate).`


    Parameters
    ----------
    sdate : date like, optional
        Calendar start date. The default is `'1980-01-01'`.
    edate : date like, optional
        Calendar end date. The default is `'2050-12-31'`.
        
    Returns
    -------
    `numpy.busdaycalendar` : NYSE business calendar.
    """
    warnings.warn('\nNYSEgen will be deprecated in future versions\n' + 
                  'instead please use calendarGen(name="NYSE", sdate=sdate, edate=edate)\n')
    
    sdate_ = np.datetime64(sdate)
    edate_ = np.datetime64(edate)
    # get the NYSE holiday list from standard pandas_market_calendars
    # ldates = list(mcal.get_calendar('NYSE').holidays().holidays)
    # append here additional holidays 
    # already added to the new version of pandas_market_calendars 
    # ldates.append(np.datetime64('2012-10-29'))
    # ldates.append(np.datetime64('2012-10-30'))
    # ldates.append(np.datetime64('2001-09-11'))
    # ldates.sort()
    # move to np.array and return the business calendar
    # hdates = np.array(ldates)
    hdates = np.array(mcal.get_calendar('NYSE').holidays().holidays)
    hdates = hdates[(hdates >= sdate_) & ( hdates <= edate_)]
    
    return np.busdaycalendar(holidays=hdates)


def calendarGen(name='NYSE', sdate='1980-01-01', edate='2050-12-31'):
    """
    Returns exchange business calendar

    Parameters
    ----------
    name : `str`, optional
        The exchange name. A valid exchange name is listed \n
            `get_calendar_names()` \n
        Default value is `'NYSE'` - i.e. New York Stock Exchange.
    sdate : date like, optional
        Calendar start date. The default is `'1980-01-01'`.
    edate : date like, optional
        Calendar end date. The default is `'2050-12-31'`.
        
    Returns
    -------
    `numpy.busdaycalendar` : NYSE business calendar.
    """
    sdate_ = np.datetime64(sdate)
    edate_ = np.datetime64(edate)
    hdates = np.array(mcal.get_calendar(name).holidays().holidays)
    hdates = hdates[(hdates >= sdate_) & ( hdates <= edate_)]

    return np.busdaycalendar(holidays=hdates)


def get_calendar_names():
    """
    Returns calendar exchange names.
    """
    return mcal.get_calendar_names()


def _set_calendar_exchange(calendar):
        if calendar is None:
            return calendarGen()
        elif isinstance(calendar, str):
            try:
                return calendarGen(calendar)
            except:
                raise Exception("Unknown exchange calendar name")
        elif isinstance(calendar, np.busdaycalendar):
            return calendar
        else:
            raise Exception("Unknown calendar type. Must be str (valid name of an exchange calendar), " +
                            "numpy.busdaycalendar, or None")