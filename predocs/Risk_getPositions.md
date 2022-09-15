#### <span style="color:green">getPositions</span>

Computes the rebalanced and delta numbers of shares for each portfolio
component.

*Call:*

```
getPositions(nshares=None, cash=0, ww=None, rtype=None, mu=None, mu0=0.,
             aversion=None, ww0=None, )
```

*Inputs:*

* `nshares` : `panda.Series`, optional;
    Initial number of shares per portfolio component.
    A missing component
    entry will be considered 0. A `None` value assumes that all
    components entries are 0. The name of the components must be
    present in the mrkdata. The default is `None`.
* `cash` : `float`, optional;
    Additional cash to be added to the capital. A
    negative entry assumes a reduction in the total capital
    available for rebalance. The total capital cannot be < 0.
    The default is `0`.
* `ww` : `panda.Series`, optional;
    External overwrite portfolio weights.
    If it not set to `None` these
    weights will overwrite the calibrated.
    The default is `None`.
* `rtype` : `str`, optional;
    Optimization type. If is not `None` it will overwrite the value
    set by the constructor. The default is `None`.
* `mu` : `float`, optional
    Targeted portfolio expected rate of return.
    Relevant only if `rtype='Risk'`
    The default is `None`.
* `mu0` : `float`, optional;
    Risk-free rate accessible to the investor.
    Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
    The default is `0`.
* `aversion` : `float`, optional;
    The value of the risk-aversion coefficient.
    Must be positive. Relevant only if `rtype='RiskAvers'`.
    The default is `None`.
* `ww0` : `list` (also `numpy.array` or `pandas.Series`), optional;
    Targeted portfolio weights
    Relevant only if `rype='InvNrisk'`.
    Its length must be equal to the number of
    symbols in market data.
    All weights must be >= 0 with sum > 0.
    If it is a list or a `numpy.array` then the weights are assumed to
    by in order of rrate.columns. If it is a `pandas.Series` then the index
    should be compatible with the `rrate.columns` or mktdata symbols
    (same symbols, not necessary in the same order).
    If it is `None` then it will be set to equal weights.
    The default is `None`.

*Returns:* `pandas.DataFrame` containing the rolling information.

Columns:

* `'old_nsh'` : initial number of shares per portfolio component as well as
additional cash position. These are present in the input.
* `'new_nsh'` : the new number of shares per component plus the  residual
cash (due to the rounding to an integer number of shares). A negative entry
means that the investor needs to add more cash in order to cover for the
number of share  roundups. It has a small value.
* `'diff_nsh'` : delta number of shares - the number of shares that needs to be
both/sold in order to rebalance the portfolio positions.
* `'weights'` : portfolio weights used for rebalance. The `'cash'` entry
is the new portfolio value (invested capital).
* `'prices'` : share prices used for rebalance evaluations.

>Note: Since the prices are closing prices, the rebalance can be executed next
business day. Additional cash slippage may occur due to share price differential
between the previous day closing and  execution time.
