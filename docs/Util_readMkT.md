

# Read historical market data <a name="TOP"></a>

There are 2 ways to get historical time series. Using:
1. [`readMkT`](#readMkT) function. It is a coinvent wrapper around `MkTreader`
class,
2. [`MkTreader`](#MkTreader) class.

In both cases the user has the abilities to:
* save the collected historical time
  series to a local file repository for later use,
* read directly from the providers,
* read and update an existing local repository.

The following market providers can be accesses:
* *yahoo* - free *as is* service,
* *eodhistoricaldata* - needs a premium account (some tests may be accomplished
  with a free key),
* *alphavantage* - needs a premium account,
* *marketstack* - needs a premium account (some tests may be accomplished
  with a free key),
* *eodhistoricaldata_yahoo* - hybrid service where the historical raw prices
are collected from *eodhistoricaldata* provider while the splits
and dividends are collected from *yahoo*. It may work with a
_eodhistoricaldata_ free key for limited testing purposes.
* *alphavantqge_yahoo* - hybrid service where the historical raw prices
are colleected from *alphavantage* provider while the splits
and dividends are collected from *yahoo*. It may work with a
_alphavantage_ free key for testing purposes.

Regardless where the historical data is collected, the
`open`, `high`, `low`, `close`, `volume` `divd` (dividends) are splits adjusted
while the close adjusted, for short `adjusted` prices, are splits and
dividends adjusted relative to the most recent date in the time series.

The following file formats are supported to save data:
* *csv* - comma separated values format
* *json* - JavaScript object notation format
* *feeder* - portable binary format for storing Arrow tables and data frames in
Python and R

In addition of the marker information, the source and acquisition date and time
are also saved.

<a name="readMkT">

## **readMkT function**
It is a convenient wrapper of `MkTreader` class. It returns an object
(`pandas.DataFrame` or `dict` of `pandas.DataFrames`) containing the
requested historical data for a set of stock symbols.

*Call:*

```
readMkT(symbol=[], sdate="2012-01-01", edate='today', calendar=None,
        output_format='frame', source=None, force=False, save=True,
        file_dir="outDir", file_format='csv', api_key=None, param=None,  
        verbose=True)
```

*Inputs:*

* `symbol` : `str, list`. <br>
A `str` containing a single stock symbol or a list of stock symbols.
The default is `[]`.
* `sdate` : `date like`. <br>
Start date of historical time series requested.
If `sdate` is not a valid exchange business day then it will
be adjusted forward to the next business day.
The default is `'2012-01-01'`.   
* `edate` : `date like`. <br>
End date of historical time-series requested. This is the date
relative to which the adjusted data is computed. If `edate` is not a
valid exchange business day then it will be adjusted backwards to the
pervious business day.
The `edate` must be greater or equal to the `sdate`.
The default is the current date, `'today'`.
* `calendar` : `numpy.busdaycalendar`. <br>
Is the exchange business calendar. If it is set to `None` then the
NY stock exchange (as returned by `azapy.NYSEgen()` function) is
assumed. The default value is `None`.
* `output_format` : `str`. <br>
The function output format. It could be
  - `'frame'` : the output format is a `pandas.DataFrame` where the
  _index_ `'date'` are the historical observation dates and the *columns* are
  '`symbol`', `'open'` , `'high'`, `'low'`, '`close`', '`volume`', '`adjusted`',
  '`divd`' and `'split'`.
  - `'dict'` : the output is a `'dict'` where the *keys* are the stock
  symbols and the *values* are `pandas.DataFrame` in the same format as
  above.

  The default is '`frame`'
* `source` : `str, dict`. <br>
If it is a `str` then it indicates a the market data provider (`'yahoo'`,
`'eodhistoricaldata'`, `'alphavantage'`, `'marketstack'`,
'`eodhistoricaldata_yahoo'`, `'alphavantage_yahoo'`). A value of `None` will
assume *yahoo* service. <br>
It can be set as a `dict` in order to include specific information for
each stock symbol. The *keys* of the dictionary are stock symbols, while then
_values_ are `dict` with additional information. The *keys* for the
information `dict` are (same as the names of the calling variables):
`source`, `force`, `save`, `file_dir`, `file_foramt`, `api_key`,
`param` and `verbose`. An absent *key* will be
filled with the value provided by the corresponding variable from the function
call. The actual set of symbols in the request is the union of
`symbol` and the *keys* of the `dict source`. <br>
Example: <br>
```
  source = {'AAPL': {'source': 'eodhistoricaldata', 'force': True},
            'SPY': {'source': 'yahoo', 'save': False}}
```
In this case there are only 2 symbol listed. They will be added (union) to
list from `symbol`. For first symbol, `'AAPL'`, the market data
provider is *eodhistoricaldata* while for `'SPY'` is *yahoo*. Also
for `'AAPL'` the flag `force` is set to `True` while for `'SPY'` the
flag `save` is set to `False`. In each case the rest of the information
is set to the values of call variables.
* `force` : Boolean flag. <br>
  - `True` : market data will be read directly from the provider.
  - `False` : it will try first to read the market data from the local
  saved file. If the saved file is not funded than it will read from
  the provider. If the saved information is too short then the missing data
  will be updated from the provider.

  The default is `False`.
* `save` : Boolean flag. <br>
  - `True` : market data will be saved in a local file.
  - `False` : suppress any saving of data.

  The default is `True`.
* `file_dir` : `str`. <br>
Local directory where market data can be read/save. If the directory does
not exist then it will be created. <br>
The default is `'outData'`.
* `api_key` : `str`. <br>
The market data provider API key where is necessary. If it is set to
`None` then the key will be read from the global environment
variables:
  - `EODHISTORICALDATA_API_KEY` - for *eodhistoricaldata.com* provider,
  - `ALPHAVANTAGE_API_KEY` - for *alphavantage.co* provider,
  - `MARKETSTACK_API_KEY` - for *marketstack.com* provider.

  The default is `None`.
* `file_format` : `str`.
The saved data files format. It could be: '`csv`', '`json`' or '`feather`'.
If it is set to `None` the `'csv'` format is assumed. <br>
> The full name of the file is: `file_dir + symbol + '.' + file_format`

  The default is `None`.
* `param` : `dict`.
Additional information to access the market data provider.
The *alphavantage* provider requires the maximum number of requests
per minute parameter. Its vale varies with the account (key) level
of access. The minimum value is 5 for a free key
and starts at 75 for a premium key. This value is stored under
`'max_req_per_min'` `dict` key. <br>
Example: `param = {'max_req_per_min': 5}` <br>
This is also the default vale for *alphavantage*, if `param` is set to
`None`. <br>
The default is `None`.
* `verbose` : Boolean flag. <br>
  - `True`: print the progress information.
  - `False`: suppress the printing of progress information.

  The default is `True`.

*Returns:*

A `pandas.DataFrame` or a dictionary of `pandas.DataFrame` according to
the value set by the `output_format` variable.

>Hints:
- use the flag combination `force=True` and `save=True` to overwrite the
an existing old data file in `file_dir`.
- the combination `force=False` and `save=True` will update an existing
old data file present in `file_dir`.

[TOP](#TOP)

### Examples:

```
import azapy as az

sdate = "2000-01-01"
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "../../MkTdata"

# simple calls
# returns a pd.DataFrame
mktdata = az.readMkT(symb, sdata=sdate, sdata=edate, file_dir=mktdir)

# returns a dict of pd.DataFrame
mktdata_dict = az.readMkT(symb, sdata=sdate, edata=edate, file_dir=mktdir,
                          output_format='dict')

# complex call
source = {'GLD': {'force': True,
                  'save': True,
                  'file_dir': '../../MkTdata_yahoo',
                  'foramt_foramt' : 'feather'
                 }
          'TLT': {'source' : 'alhavantage',
                  'force': False,
                  'save': True,
                  'file_dir': '../../MkTdata_av',
                  'foramt_foramt': 'json'
                  'param': {'max_req_per_min': 75}
                 }
          'XLV': {'source': 'eofhistoricaldata'}
         }

symb = ['VGT', 'PSJ']

mktdata = az.readMkT(symb, sdata=sdate, sdata=edate, source=source,
                     file_dir=mktdir, file_format='csv')
```
[TOP](#TOP)

<a name="MkTreader">

## **class MkTreader**

Process historical market data requests.

**Methods:**

* [<span style="color:green">get</span>](#get)
* [<span style="color:green">get_request_status</span>](#get_request_status)
* [<span style="color:green">get_error_log</span>](#get_error_log)

### Constructor
```
  MkTreader()
```
There are no variables to be passed to the constructor.

[TOP](#TOP)

### Methods:

<a name="get">

#### <span style="color:green">get</span>

It is the main function that executes the data request.

*Call:*

```
get(symbol=[], sdate="2012-01-01", edate='today', calendar=None,
    output_format='frame', source=None, force=False, save=True,
    file_dir="outDir", file_format='csv', api_key=None, param=None,  
    verbose=True)
```

*Inputs:*

The inputs are the same as for [`readMkT` function](#readMkT).

[TOP](#TOP)

<a name="get_request_status">

#### <span style="color:green">get_request_status</span>

Returns the request status per symbol.

*Call:*

```
get_request_status()
```

*Inputs:*

There are no input variables.

*Returns:*

The output is a `pandas.DataFrame` contained, per symbol, the actual request
parameters (except the API keys if they were read from the global environment
variables), as well as the `'error'` status. A value of `'No'` it means that
there are some missing observation dates.
The list of missing observation dates can be further explore by calling the
member function `get_error_log()`.
A `'Yes'` means that everything is OK.
A value of `None` occurs if the request was not completed (other warning
messages may clarify these situations).

[TOP](#TOP)

<a name="get_error_log">

#### <span style="color:green">get_error_log</span>

Returns the lists of missing observations dates per symbol.

*Call*

```
get_error_log()
```

*Inputs:*

There are no input variables.

*Returns:*

The return is a `dict` where the _keys_ are the names of the symbols with
missing observation dates. If the `dict` is empty, then there are no missing
observations for any symbol. The _values_ of the `dict` are also `dict` with
possible _keys_:
- `'front'` : the list of missing dates between the most recent observation
date in the time series and the request `edate` (adjusted backward to a
business day),
- `'back'` : the list of missing dates between the request  
`sdate` (adjusted forward to a business day) and the oldest observation
date in the time series,
- `'mid'` : the list of missing dates between the oldest and the most recent observation dates of the time series.

Only the keys with non-empty lists of missing dates will be present.

[TOP](#TOP)

### Examples:

```
import azapy as az

sdate = "2000-01-01"
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "../../MkTdata"

mkt = az.MkTreader()
hdata = mkt.get(symb, sdata=sdate, sdata=edate, file_dir=mktdir)

# request status
req_status = mkt.get_request_status()

# missing observation dates
error_date = mkt.get_error_log()
```

[TOP](#TOP)
