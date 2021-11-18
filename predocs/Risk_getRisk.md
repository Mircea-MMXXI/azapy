
<a name="getRisk"></a>

#### <span style="color:green">getRsik</span>

Computes the risk of a portfolio defined by a set of weights.

*Call:*
```
getRisk(ww, rrate=None)
```

*Inputs:*

* `ww` : List like of portfolio weights. Its length must be equal to the
number of symbols in `rrate` (`mktdata`). All weights must by `>0`. If it
is a `list` or a `np.array` then the weights are assumed to be in the order
of `rrate.columns`. If it is a `pd.Series` the index should be compatible
with the `rrate.columns` or `mktdata` symbols (not necessary in the same
order).
* `'rrate'` : `pd.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed by the constructor from `mktdata`. The default is `None`.

*Returns:* The value of the risk measure.

[TOP](#TOP)

---
