#### <span style="color:green">getRisk</span>

Computes the risk of a portfolio defined by a set of weights.

*Call:*
```
getRisk(ww, rrate=None)
```

*Inputs:*

* `ww` : `list-like` of portfolio weights. Its length must be equal to the
number of symbols in `rrate` (`mktdata`). All weights must by $\ge 0$ and
their sum equal to $1$. If it
is a `list` or a `numpy.array` then the weights are assumed to be in the order
of `rrate.columns`. If it is a `pandas.Series` the index should be compatible
with the `rrate.columns` or `mktdata` symbols (not necessary in the same
order).
* `'rrate'` : `pandas.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed by the constructor from `mktdata`. The default is `None`.

*Returns:* The value of the risk measure.

Note: It will set the following class members:
* _risk_
* _primary_risk_comp_
* _secondary_risk_comp_
* _sharpe_
* _RR_
