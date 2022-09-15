#### <span style="color:green">set_rrate</span>

Sets portfolio components historical rates of returns.
It will overwrite the value computed by the constructor from `mktdata`.

*Call:*

```
set_rrate(rrate)
```

*Inputs:*

* `rrate` : `pandas.DataFrame`;
Portfolio components historical rates of returns. The index name is `'date'` and
columns are `symbol1`, `symbol2`, etc.


*Returns:* `None`
