
<a name="getRisk">

#### <span style="color:green">getRsik</span>

Computes the risk of a portfolio

*Call:*
```
getRisk(ww, rrate=None)
```

*Input:*

* `ww` : List like of portfolio weights. Its length must be equal to the
number of symbols in `rrate` (`mktdata`). All weights must by $>0$. If it
is a `list` or a `np.array` then the weights are assumed to by in order
of `rrate.columns`. It is a `pd.Series` the index should be compatible
with the `rrate.columns` or `mktdata` symbols (not necessary in the same
order).
* `'rrate'` : `pd.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed in the constructor from `mktdata`. The default is `None`.

*Returns:* The value of the risk measure.

[TOP](#TOP)

---
