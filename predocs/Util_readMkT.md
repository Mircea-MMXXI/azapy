
# readMkT

## Collects historical market data from *alphavantage* provider.

A valid key is needed to access *alphavantage* servers
(*see alphavantage.co for more information on how to obtain a key*).
The key value can be stored in the environment variable
'ALPHAVANTAGE_API_KEY' or it can be passed directly to the `readMkT` function.

The function offers the facility to save the data on the local drive for future
use. It is intended to speed data access in working sessions requiring multiple
readings of the same data and to reduce the number of external
servers hits.

### Call:

```
readMkT(symbols,
        dstart=np.datetime64('2012-01-01'), dend=np.datetime64('today'),
        force=False, verbose=True, dir='outData', save=True,
        api_key=None, maxsymbs=5, adj_split=True, out_dict=False)
```

### Inputs:
* `symbols` :
    List of stock symbols.
* `dstart` : Start date of historical time-series. The default is
`np.datetime64('2012-01-01')`.
* `dend` : End date of historical time-series. The default is the current
date, `np.datetime64('today')`.
* `force` : Boolean flag.
    - `True` : market data will be read from *alphavantage*.
    - `False` : it attempts first to read the market data from the local
    directory indicated by the variable `dir`. If it
    fails, then it will read from *alphavantage*.

    The default is `False`.
* `verbose` : Boolean flag.
    - `True`: print the progress information.
    - `False`: suppress the printing of progress information.

    The default is `True`.
* `dir` :
    Local directory where market data can be read/save.
    The default is 'outData'.
* `save` : Boolean flag.
    - `True` : market data will be saved to `dir`.
    - `False` : suppress any saving of data.

    The default is `True`.
* `api_key` :
    Valid *alphavantage* key. If is set to `None` then the key will be read
    from environment variable 'ALPHAVANTAGE_API_KEY'.
    The default is `None`.
* `maxsymbs` :
    Maximum number of symbols that can be read per minute. Should match
    *alphavantage* account limits. After the `maxsymbs` is reached the program
    will sleep for 1 min before resuming reading more data.
    The default is 5.
* `adj_split` : Boolean flag.
    - `True` : adjust all fields for split events.
    - `False` : no adjustment is made.

    The default is `True`.
* `out_dict` : Boolean flag.
    Choose the output format:

    - `True` : `pd.DataFrame` with columns: "symbol", "date", "open",
    "high", "low", "close", "volume", "adjusted", "divd", "split".

    - `False` : a dictionary where the keys are the symbols and the items are
    the `pd.DataFrame` with columns: "date", "open", "high",
    "low", "close", "volume", "adjusted", "divd", "split".

    The default is `False`.

### Returns:

A `pd.DataFrame` or a dictionary of `pd.DataFrame` according to the value set of
`out_dict`.


>Hint: use the flag combination `force=True` and `save=True` to overwrite the
existing old data from `dir`.

### Examples:

```
import pandas as pd

import azapy as az

#==============================================================================
# Collect market data
sdate = pd.to_datetime("2000-01-01")
edate = pd.to_datetime('today')
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "../../MkTdata"

# force=True read directly from alphavantage
# force=False read first from local directory, if data does not exists,
#             read from alphavantage

# returns a pd.DataFrame
mktdata = az.readMkT(symb, dstart = sdate, dend = edate,
                     dir=mktdir, force=False)

# returns a dict of pd.DataFrame
mktdata_dict = az.readMkT(symb, dstart=sdate, dend=edate, force=False,
                          dir=mktdir, out_dict=True)
```
