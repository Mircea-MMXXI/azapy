import pandas as pd
import numpy as np
from collections import defaultdict 

def prep_uw_(rdata):
    """
    To be used internaly
    Computes the underwater vector
    
    Parameters
    ----------
    rdata : pd.Series
        Series of price indexed by dates

    Returns
    -------
    TYPE: pd.Sereis
        The underwater vector
    """
    return rdata / rdata.cummax() - 1


def max_drawdown_(uw):
    """
    To be used internally
    Compute the maximum drawdown for an underwater vector

    Parameters
    ----------
    uw : TYPE dp.Series
        under-water vector procuded by prep_uw_ function

    Returns
    -------
    draw_val : TYPE float
        the value of the drawdown
    draw_min : TYPE pd.Timestamp
        date of the maximum drawdown
    draw_start : TYPE pd.Timestamp
        date of the drawdown start
    draw_end : TYPE pd.Timestamp
        date of the drawdown recovery 
        if in np.nan then the drawdown is in progress
    """
    draw_min = uw.idxmin()
    draw_val = uw[draw_min]
    draw_start = uw[:draw_min][uw[:draw_min] == 0].index[-1]
    try:
        draw_end = uw[draw_min:][uw[draw_min:] == 0].index[0]
    except IndexError:
        draw_end = np.nan  # drawdown not recovered
        
    return draw_val, draw_min, draw_start, draw_end

def max_drawdown(rprice, col=np.nan):
    """
    Compute the max_drowdown for a price time-series

    Parameters
    ----------
    rprice : TYPE pd.Series or pd.DataFram
        time-sereis of prices as a pd.Series or as column in a dp.DataFrame
    col : TYPE, string
        column name if rprice is a DataFrame. If is set to np.nan then rprice
        is assumed to be a Series. The default is np.nan.

    Returns
    -------
    val : TYPE float
        The value of the drawdown.
    i_min : TYPE pd.Timestamp
        The maximum drawdown date.
    i_start : TYPE pd.Timestamp
        Date when the drawdown had started.
    i_end : TYPE pd.Timestamp
        Date of the drawrown recovery. A value of np.nan indicates that the
        drawdwon is in progress.
    """
    rdata = rprice if pd.isna(col) else rprice[col]

    val, i_min, i_start, i_end = max_drawdown_(prep_uw_(rdata))
    
    i_min = i_min.strftime('%Y-%m-%d')
    i_start = i_start.strftime('%Y-%m-%d')
    if not pd.isna(i_end): 
        i_end = i_end.strftime('%Y-%m-%d')
        
    return val, i_min, i_start, i_end

def drawdown(rprice, col=np.nan, top=10):
    """
    Computes the largest drawdowns for a price time-series

    Parameters
    ----------
    rprice : TYPE pd.Series or pd.DataFrame
        time-sereis of prices as a pd.Series or as column in a dp.DataFrame
    col : TYPE, string
        Name of the column of price if rprice is a pd.DataFrame. If its value 
        is set to np.nan then rprice is assumed to be a pd.Series.
        The default is np.nan.
    top : TYPE, int
        Maximum number of the largest drawdown to be computed. 
        The default is 10.

    Returns
    -------
    TYPE pd.DataFrama
        Tabel containing the drawdowns ordered from the largest to smallest.
        Tabel columns are:
                'DD': (folat) drawdown max value
                'Date': (pd.Timestamp) drawdown max value date
                'Start': (pd.Timestamp) drawdown start date
                'End': (pd.Timestampo) drawdown recovery date
        The number of rows is <= top (smaller then top if max number of 
                                      drawdons is smaller than top)
    """
    rdata = rprice if pd.isna(col) else rprice[col]
 
    uw = prep_uw_(rdata)
    dd = defaultdict(lambda: [])
    
    for _ in range(top):
        val, i_min, i_start, i_end = max_drawdown_(uw)
        
        dd['DD'].append(val)
        dd['Date'].append(i_min.strftime('%Y-%m-%d'))
        dd['Start'].append(i_start.strftime('%Y-%m-%d'))
        if not pd.isna(i_end): 
            dd['End'].append(i_end.strftime('%Y-%m-%d'))
        else:
            dd['End'].append(np.nan)
            i_end = uw.index[-1]
    
        uw.drop(index=uw[i_start:i_end].index[:], inplace=True)
        if (len(uw) == 0) or (uw.min() == 0): break

    return pd.DataFrame(dd, 
                        index=pd.Index(range(1, len(dd['DD']) + 1), name='No'))

# Examples
# from readMkTData import readMkT

# symb = ['AAPL']
# rprice = readMkT(symb, force=False, adj_split=False, out_dict=False)

# max_drawdown(rprice, col='adjusted')
# max_drawdown(rprice['adjusted'])
# %timeit max_drawdown(rprice['adjusted'])
# drawdown(rprice, col='adjusted', top=5)
# drawdown(rprice['adjusted'], top=15)

# rprice['adjusted'].plot()
# rprice.loc[rprice.index >= pd.to_datetime("2021-01-01"), 'adjusted'].plot()
# (rprice.loc[pd.to_datetime("2020-03-23"),'adjusted'] / 
#   rprice.loc[pd.to_datetime("2020-02-12"),'adjusted'] - 1)

# rprice['cash'] = 1
# max_drawdown(rprice['cash'])
# drawdown(rprice['cash'])

# rprice['rinfl'] = 0.03 / 254
# rprice['infl'] = (rprice['rinfl'] + 1).cumprod()
# max_drawdown(rprice['infl'])
# drawdown(rprice['infl'])

# rprice['rdinfl'] = -0.03 / 254
# rprice['dinfl'] = (rprice['rdinfl'] + 1).cumprod()
# rprice.dinfl[-1] / rprice.dinfl[0] - 1
# max_drawdown(rprice['dinfl'])
# drawdown(rprice['dinfl'])


