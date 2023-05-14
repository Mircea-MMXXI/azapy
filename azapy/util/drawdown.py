""""""  # <- add this
"""
Contains:
    - max_drawdown : returns maximum drawdown
    - drawdown : returns a list of drawdowns
"""
import pandas as pd
import numpy as np
from collections import defaultdict

def _prep_uw(rdata):
    """
    Computes the underwater vector
    """
    return rdata / rdata.cummax() - 1


def _max_drawdown(uw):
    """
    Computes the maximum drawdown for an underwater vector
    """
    draw_min = uw.idxmin()
    draw_val = uw[draw_min]
    draw_start = uw[:draw_min][uw[:draw_min] == 0].index[-1]
    try:
        draw_end = uw[draw_min:][uw[draw_min:] == 0].index[0]
    except IndexError:
        draw_end = np.nan  # drawdown not recovered

    return draw_val, draw_min, draw_start, draw_end

def max_drawdown(mktdata, col=None):
    """
    Computes the maximum drawdown for a time-series of prices.

    Parameters
    ----------
    mktdata : `pandas.Series` or `pandas.DataFram`
        Time-series of prices as a `pandas.Series` or as column in a 
        pandas.DataFrame`
    col : `str`, optional
        Column name if `mktdata` is a `pandas.DataFrame`. If is set to `None`, 
        then `mktdata` is assumed to be a `pandas.Series`. 
        The default is `None`.

    Returns
    -------
    (`float`, `pandas.Timestamp`, `pandas.Timestamp`, `pandas.Timestamp`) : Tuple 
        - value of the drawdown, 
        - maximum drawdown date, 
        - drawdawn start date,
        - drawdown end date. It is set to `nan` if the drawdown is still 
          in progress.
    """
    rdata = mktdata if col is None else mktdata[col]

    val, i_min, i_start, i_end = _max_drawdown(_prep_uw(rdata))

    i_min = i_min.strftime('%Y-%m-%d')
    i_start = i_start.strftime('%Y-%m-%d')
    if not pd.isna(i_end):
        i_end = i_end.strftime('%Y-%m-%d')

    return val, i_min, i_start, i_end

def drawdown(mktdata, col=None, top=10):
    """
    Computes the largest drawdowns for time-series of prices.

    Parameters
    ----------
    mktdata : `pandas.Series` or `pandas.DataFrame`
        Historical daily time-series of prices as a `pandas.Series` 
        or as a `pandas.DataFrame`
    col : `str`, optional
        The column name if mktdata is a `pandas.DataFrame`. 
        If it is set to `None`, then `mktdata` is assumed to be a 
        `pandas.Series`.
        The default is `None`.
    `top` : `int`, optional
        Maximum number of largest drawdowns to be reported.
        The default is `10`.

    Returns
    -------
    `pandas.DataFrama` : Table containing the drawdowns ordered from the 
    largest to smallest.
        Table columns are:

            - 'DD': (float) drawdown max value
            - 'Date': (`pandas.Timestamp`) drawdown max value date
            - 'Start': (`pandas.Timestamp`) drawdown start date
            - 'End': (`pandas.Timestamp`) drawdown recovery date
            
        The number of rows is <= top
    """
    rdata = mktdata if pd.isna(col) else mktdata[col]

    uw = _prep_uw(rdata)
    dd = defaultdict(lambda: [])

    for _ in range(top):
        val, i_min, i_start, i_end = _max_drawdown(uw)

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
