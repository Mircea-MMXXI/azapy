
<a name="set_model">

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.


*Call:*

```
set_model(ww=None)
```

*Input:*

* `ww` :
List (alos `np.array` or `pd.Series`) of weights. If it is a `pd.Series`
the index should match the basket `symb`.
Otherwise the weights are considered in the  `symb`
order. If it is set to `None` than `ww` will be set to equal weights.
The default is `None`.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---
