There are 2 ways to get historical time series:
1. [**readMkT**](azapy.MkT.readMkT.readMkT) function.
   It is a coinvent wrapper around `MkTreader` class,
2. [**MkTreader**](azapy.MkT.MkTreader.MkTreader) class.
   Provides additional facilities for tracking missing data and errors.

In both cases the user can:
* save the collected historical time
  series to a local file repository for later use,
* read directly from the providers,
* read and update an existing local repository.

The following market providers can be accesses:
* <span style="color:green">**yahoo**</span> - free *as is* service,
* <span style="color:green">**eodhistoricaldata**</span> - needs a premium
account from *eodhistoricaldata.com*
(some tests may be accomplished with a free key),
* <span style="color:green">**alphavantage**</span> - needs a premium account
from *alphavantage.co*,
* <span style="color:green">**marketstack**</span> - needs a premium account
from *marketstack.com*
(some tests may be accomplished with a free key),
* <span style="color:green">**eodhistoricaldata_yahoo**</span> - hybrid service
where the historical raw prices
are collected from *eodhistoricaldata* provider while the splits
and dividends are collected from *yahoo*. It may work with a
_eodhistoricaldata_ free key for limited testing purposes.
* <span style="color:green">**alphavantage_yahoo**</span> - hybrid service
where the historical raw prices
are collected from *alphavantage* provider while the splits
and dividends are collected from *yahoo*. It may work with a
_alphavantage_ free key for testing purposes.

Regardless of where the historical data is collected, the
`open`, `high`, `low`, `close`, `volume`, `divd` (dividends) are splits adjusted
while the close adjusted, for short `adjusted` prices, are splits and
dividends adjusted relative to the most recent date in the returned time series.

The following file formats are supported to save data:
* *csv* - comma separated values format
* *json* - JavaScript object notation format
* *feeder* - portable binary format for storing Arrow tables and data frames in
Python and R

In addition of the marker information, the source and acquisition date and time
are also saved.
