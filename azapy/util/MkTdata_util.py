"""
Contains:
    
    - add_cash_security : and a cash like security
    - update_all_MkTData : update mkt data saved in a directory
"""

import pandas as pd
import os
import warnings
from azapy.MkT.MkTreader import MkTreader

def add_cash_security(data, name='_CASH_', value=1.):
    '''
    Add a cash like security to the MkT data.

    Parameters
    ----------
    data : pandas.DataFrame or dict
        MkT data to add a cash like positions.
    name : str, optional
        The symbol of the cash like security. Must be different than any 
        symbol in data. Note that CASH is a valid stock symbol.
        The default is '_CASH_'.
    value : float, optional
        Nominal value of the cash like security. It is 
        constant over time (no dividends and no splits).
        The default is 1.

    Returns
    -------
    pandas.DataFrame (if data was a pandas.DataFame)
        The new MkT data updated with the cash like security.
    dict (if data was a dict)
        Append to initial MkT data dict a cash like security.
    '''
    if isinstance(data, pd.core.frame.DataFrame):
        cash = pd.DataFrame(value, index=data.index.unique(), 
                            columns=data.columns)
        cash['symbol'] = name
        cash['divd'] = 0.
        return pd.concat([data, cash])
    
    elif isinstance(data, dict):
        ld = list(data.values())[0]
        cash = pd.DataFrame(value, index=ld.index, columns=ld.columns)
        cash['symbol'] = name
        cash['divd'] = 0.
        data[name] = cash
        return data
    
    else:
        raise TypeError(f"data: wrong type {type(data)} "
                        + "must be pandas.DataFrame or dict")
        
    
def update_all_MkTData(mktdir, source=None, api_key=None, param=None,
                       except_file=[], verbose=True):
    '''
    Updates all mkt data saved in a directory.

    Parameters
    ----------
    mktdir : str
        Mkt data directory.
    source : str, optional
        Mkt data provider. 
        For more details see the `azapy.MkTreader.get` function doc. 
        The default is `None`.
    api_key : str, optional
        Mkt data provider API key.
        For more details see the `azapy.MkTreader.get` function doc. 
        The default is `None`
    param : dict, optional
        Additional parameters required by mkt data provider.
        For more details see the `azapy.MkTreader.get` function doc. 
        The default is `None`.
    except_file : list, optional
        List of symbols to be omitted from the update. The default is [].
    verbose : Boolean, optional
    
        - `True` will print a progress report,
        - `False` suppress any printing to the terminal.
        
        The default is `True`.

    Returns
    -------
    Error code (int):
        - 200 : successful, everything updated
        - 201 : some (or all) were not completly updated 
        - 101 : the `mktdir` does not exists
        - 102 : unsupported mkt data sources
            
    Note that files with unsupported extensions 
    (see `azapy.MkTreader.get` function)
    are silently omitted from the update.
    '''
    
    if not os.path.isdir(mktdir):
        warnings.warn(f"Warning: doesn\'t exist dir: {mktdir}")
        return 101
    if source is None:
        source = 'yahoo'
    elif not source in ['yahoo', 
                        'eodhistoricaldata', 'eodhistoricaldata_yahoo',
                        'alphavantage', 'alphavantage_yahoo',
                        'marketstack']:
        warnings.warn(f"Unsupported mkt data source {source}")
        return 102
    if except_file is None:
        except_file = []
    
    sce = {}
    for file in os.listdir(mktdir):
        if verbose:
            print(file)
        fname, fext = os.path.splitext(file)
        if fext[1:] not in {'csv', 'json', 'feather'}:
            if verbose:
                print(f"file {file} omitted from update")
            continue
        if fname in except_file:
            if verbose:
                print(f"file {file} excepted from update")
            continue
        ld = {}
        ld['source'] = source
        ld['file_format'] = fext[1:]
        sce[fname] = ld

    mkt = MkTreader()
    mkt.get(source=sce, file_dir=mktdir, api_key=api_key, 
            param=param, verbose=verbose)
    
    elo = mkt.get_error_log()
    if verbose:
        sta = mkt.get_request_status()
        print("\nThe following symbols were updated:")
        for symb in sta.columns:
            print(f"symb {symb} : errors {sta.loc['error', symb]}")
        if len(elo) != 0:
            print(f"errors: {elo}")
          
    if len(elo) == 0: 
        return 200
    return 201