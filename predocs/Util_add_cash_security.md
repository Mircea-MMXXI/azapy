
# add_cash_security

## Add a cash like security to an existing MkT data object.

### Call:

```
  add_cash_security(data, name='_CASH_', value=1.)
```

### Inputs:

* `data` : pandas.DataFrame or dict <br>
A valid MkT data object as returned by the function `azapy.readMkT`.
* `name` : str <br>
The name of the cash like security. It must be different than the other
symbols in the `data`. Note that 'CASH' is a valid stock symbol.
The default is `'_CASH_'`.
* `value` : float <br>
The nominal value of the cash like security in units of currency. The
default is 1.

### Returns:
If the MkT data input, `data`, is a
* `pandas.DataFrame` then the
  function returns a new `pandas.DataFame` containing the
  original MkT data and the cash like securit,
* `dict` then the cash like security will be appended to the initial `data`
`dict`.

### [Examples:](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/util/add_cash_security_example.py)

```
import azapy as az

mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# mkt data as a DataFrame
mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

# new MkT DataFrame with a cash like security
mktdata_c = az.add_cash_security(mktdata)

# mkt data as a dict
mktdata_d = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir,
                     output_format='dict')

# append a cash like security to mktdata_d dict
mktdata_d_c = az.add_cash_security(mktdata_d)
```
