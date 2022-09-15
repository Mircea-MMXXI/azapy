
<a name="set_model"></a>

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
set_model(ww=None)
```

*Inputs:*

* `ww` :
List like positive weights, `len(ww)` must be equal to
`len(symb)`. If `ww` is a `pd.Series`
the index should match the portfolio symbols, `symb`
Otherwise the weights are considered in the  `symb`
order. If it is `None` than `ww` will be set to equal weights,
`ww = [1 / len(symb)] * len(symb)`.
The default is `None`.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---
