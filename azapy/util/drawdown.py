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
    Computes the max_drowdown for a price time-series

    Parameters
    ----------
    mktdata : pd.Series or pd.DataFram
        time-series of prices as a pd.Series or as column in a dp.DataFrame
    col : string, string
        column name if mktdata is a DataFrame. If is set to None then mktdata
        is assumed to be a Series. The default is None.

    Returns
    -------
    float
        The value of the drawdown.

            i_min : pd.Timestamp
                The maximum drawdown date.
            i_start : pd.Timestamp
                Date when the drawdown had started.
            i_end : pd.Timestamp
                Date of the drawdown recovery. A value of `nan` indicates that the
                drawdown is in progress.
    """
    rdata = mktdata if  col is None else mktdata[col]

    val, i_min, i_start, i_end = _max_drawdown(_prep_uw(rdata))

    i_min = i_min.strftime('%Y-%m-%d')
    i_start = i_start.strftime('%Y-%m-%d')
    if not pd.isna(i_end):
        i_end = i_end.strftime('%Y-%m-%d')

    return val, i_min, i_start, i_end

def drawdown(mktdata, col=None, top=10):
    """
    Computes the largest drawdowns for a price time-series

    Parameters
    ----------
    mktdata : pd.Series or pd.DataFrame
        time-series of prices as a pd.Series or as column in a dp.DataFrame
    col : string, optional
        Name of the column of price if mktdata is a pd.DataFrame. If its value
        is set to None then mktdata is assumed to be a pd.Series.
        The default is None.
    top : int, optional
        Maximum number of the largest drawdown to be computed.
        The default is 10.

    Returns
    -------
    pd.DataFrama
        Table containing the drawdowns ordered from the largest to smallest.
        Table columns are:

            - 'DD': (float) drawdown max value
            - 'Date': (pd.Timestamp) drawdown max value date
            - 'Start': (pd.Timestamp) drawdown start date
            - 'End': (pd.Timestamp) drawdown recovery date
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
