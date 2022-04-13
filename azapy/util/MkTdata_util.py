import pandas as pd

def MkTdata_add_cash(data, name='_CASH_', value=1.):
    '''
    Add a cash investment to the MkT data

    Parameters
    ----------
    data : pandas.DataFrame
        MkT data to add a cash like positions.
    name : str, optional
        The symbol of the cash like security. Must be different than any 
        symbol in data. Note that CASH is a valid stock symbol.
        The default is '_CASH_'.
    value : float, optional
        Nominal value of the cash like security. It is considered to be 
        constant over time (no dividents and no splits).
        The default is 1.

    Returns
    -------
    pd.DataFrame
        The new MkT data updated with the cash like security.
    '''
    cash = pd.DataFrame(value, index=data.index.unique(), columns=data.columns)
    cash['symbol'] = name
    cash['divd'] = 0.
  
    return pd.concat([data, cash])
