#### <span style="color:green">getWeights</span>

Computes the optimal portfolio weights.

*Call:*

```
getWeights(rtype=None, mu=None, d=1, mu0=0., aversion=None, ww0=None,
           rrate=None )
```

*Inputs:*

* `rtype` : str, optional;
    Optimization type. If is not `None` it will overwrite the value
    set by the constructor. The default is `None`.
* `mu` : float, optional.
    Targeted portfolio expected rate of return.
    Relevant only if `rtype='Risk'`
    The default is `None`.
* `d` : int, optional;
    Frontier type. Active only if `rtype='Risk'`. A value of `1` will
    trigger the evaluation of optimal portfolio along the efficient
    frontier. Otherwise, it will find the portfolio with the lowest
    rate of return along the inefficient portfolio frontier.
    The default is `1`.
* `mu0` : float, optional;
    Risk-free rate accessible to the investor.
    Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
    The default is `0`.
* `aversion` : float, optional;
    The value of the risk-aversion coefficient.
    Must be positive. Relevant only if `rtype='RiskAvers'`.
    The default is `None`.
* `ww0` : list (also `numpy.array` or `pandas.Series`), optional;
    Targeted portfolio weights.
    Relevant only if `rype='InvNrisk'`.
    Its length must be equal to the number of
    symbols in rrate (mktdata).
    All weights must be >= 0 with sum > 0.
    If it is a list or a `numpy.array` then the weights are assumed to
    by in order of `rrate.columns`. If it is a `pandas.Series` then the index
    should be compatible with the `rrate.columns` or mktdata symbols
    (same symbols, not necessary in the same order).
    If it is `None` then it will be set to equal weights.
    The default is `None`.
* `rrate` : `pandas.DataFrame`, optional;
    The portfolio components historical rates of returns.
    If it is not `None`, it will overwrite the rrate computed in the
    constructor from mktdata. The default is `None`.

*Returns:* `pandas.Series` containing the portfolio weights.

Note: It will set the following class members:
* _risk_
* _primary_risk_comp_
* _secondary_risk_comp_
* _sharpe_
* _RR_
