"""
Contains two functions

- NYSEgen : NYSE business calendar
- readMkT : read historic market data from alphavantge
"""

import pandas as pd
import pandas_datareader as web
from datetime import datetime, date
import time
import os
import pathlib
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

def readMkT(symbols, dstart=datetime(2012, 1, 1), dend=date.today(), 
            force=False, verbose=True, dir='outData', save=True,
            api_key=None, maxsymbs=5, adj_split=True, out_dict=False):
    """
    Read historic market data from alphavantage servers or from the local disk.
    Historic market data is collected between ``dstart`` and ``dend`` dates.

    Parameters
    ----------
    symbols : list
        List of stock symbols.
    dstart : datetime, optional
        Start date. The default is np.datetime64('2050-12-31').
    dend : datetime, optional
        End date. The default is date.today().
    force : boolean, optional
        - True: market data will be read from alphavantage.
        - False: it attempts to read the market data from the disk. If it
        fails then it will read from alphavantage. 
        
        The default is False.
    verbose : boolean, optional
        - True: print the progress information.
        - False: suppress the printing of progress information.
        
        The default is True.
    dir : string, optional
        Local directory where market data can be read/save. 
        The default is 'outData'.
    save : boolean, optional
        - True: market data will be saved to dir. 
        - False: suppress any saving of data. 
        
        The default is True.
    api_key : string, optional
        Valid alphavantage key. If is set to ``None`` then the key will be read
        from environment variable 'ALPHAVANTAGE_API_KEY'.
        The default is ``None``.
    maxsymbs : int, optional
        Maximum number of symbol that can be read per minute. Should match 
        alphavantage account limits. After the maxsymbs is reached the program 
        will sleep for 1 min before resuming reading more data. 
        The default is 5.
    adj_split : boolean, optional
        - True: adjust all fields for split events. 
        - False: no adjustment is made. 
        
        The default is True.
    out_dict : boolean, optional
        Choose the output format:
            
        - True: a pandas.DataFrame with columns: "symbol", "date", "open",
        "high", "low", "close", "volume", "adjusted", "divd", "split".
        
        - False: a dictoanry where the keys are the symbols and the items are
        the pandas.DataFrame with columns: "date", "open", "high",
        "low", "close", "volume", "adjusted", "divd", "split".
        
        The default is False.

    Returns
    -------
    pandas.DataFrame or a dictionary of pandas.DataFrame, according with the 
    value of ``out_dict``.
    """
    if save: pathlib.Path(dir).mkdir(parents=True, exist_ok=True) 
    rkod = 1
    rprice = {}
    for symbol in symbols:
        # Check if the symbol is already saved
        # if so then read it from the disk and go to the next
        filename = dir + "/" + symbol + ".csv"
        if not force and os.path.isfile(filename):
            sprice = pd.read_csv(filename, index_col='date')
            if verbose: print("Read " + symbol + " form " + dir)
            rprice[symbol] = sprice
            continue
        
        # Get the symbol from the web
        
        # If this is maxsymbs take 1min sleep
        if rkod > maxsymbs:
            rkod = 1
            print("Sleep 1 min.")
            time.sleep(60)
        
        # Read from web
        if verbose: 
            print("Read " + symbol + " from web")
        
        if api_key is None:
            api_key = os.getenv('ALPHAVANTAGE_API_KEY')
            
        rkod += 1
        try:
            sprice = web.DataReader(symbol, 
                                    data_source="av-daily-adjusted",
                                    start=dstart, 
                                    end=dend, 
                                    api_key=api_key)
        except:
            print("Warning: Cannot read " + symbol)
            continue
        
        # rename columns
        sprice.columns = ['open', 'high', 'low', 'close', 
                          'adjusted', 'volume', 'divd', 'split']
        sprice.index.name = 'date'
        sprice['symbol'] = symbol
        sprice = sprice[['symbol', 'open', 'high', 'low', 'close', 
                         'volume', 'adjusted', 'divd', 'split']]
        
        # Save on the disk
        if save:
            try:
                sprice.to_csv(filename, index=True)
                print( "Saved at " + filename)
            except:
                print("Cannot save " + filename)
            
        rprice[symbol] = sprice
        
    # Adjust for split
    if adj_split:
        for k, v in rprice.items():
            rprice[k] = _adjustSplit(v)
            
    # Return formatted output
    for k in rprice.keys():
        rprice[k].index = pd.to_datetime(rprice[k].index)
        
    if out_dict: 
        return rprice
    
    return pd.concat(rprice.values())
        

def _adjustSplit(data):
    spl = data.split.sort_index(ascending=False).cumprod().sort_index()
    data.open /= spl
    data.close /= spl
    data.high /= spl
    data.low /= spl
    data.volume *= spl
    data.divd /= spl
    data.split = 1
    
    return data    

    