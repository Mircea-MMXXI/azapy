
# update_all_mkTData <a name="TOP"></a>

## Update all the MkT data saved in a directory

### Call:

```
update_all_MkTData(mktdir, source=None, api_key=None, param=None,
                   except_file=[], verbose=True)
```

### Inputs:
* `mktdir` : str, <br>
Mkr data directory.
* `source` : str, <br>
Mkt data provider.
For more details see the `azapy.readMkT` function doc.
The default is `None`.
* `api_key` : str, <br>
Mkt data provider API key.
For more details see the `azapy.readMkT` function doc.
The default is `None`
* `param` : dict, <br>
Additional parameters required by mkt data provider.
For more details see the `azapy.readMkT` function doc.
The default is `None`.
* `except_file` : list, <br>
List of symbols to be omitted from the update. The default is [].
* `verbose` : Boolean, <br>

  -`True` will print a progress report,
  - `False` suppress any printing to the terminal.

The default is `True`.

### Returns:
Error code (int):

    - 200 : successful, everything updated
    - 201 : some (or all) were not completely updated
    - 101 : the `mktdir` does not exists
    - 102 : unsupported mkt data sources

> Note that files with unsupported extensions (see `azapy.readMkT` function)
are silently omitted from the update.

### [Examples:](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/util/update_all_MkTData_example.py)

```
import azapy as az

mktdir = "../../MkTdata_test"

ercode = az.update_all_MkTData(mktdir)

# error code 200 is OK (see docs)
print(f"error code: {ercode}")


```
[TOP](#TOP)
