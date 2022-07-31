#### <span style="color:green">set_rrate</span>

Sets portfolio components historical rates of returns.
It will overwrite the value computed by the constructor from `mktdata`.

*Call:*

```
set_rrate(rrate)
```

*Inputs:*

* `rrate` : `pandas.DataFrame`,
portfolio components historical rates of returns, where the
columns are `'date'`, `symbol1`, `symbol2`, etc.


*Returns:* `None`
