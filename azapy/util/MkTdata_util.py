import pandas as pd

def add_cash_security(data, name='_CASH_', value=1.):
    '''
    Add a cash like security to the MkT data

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
    
