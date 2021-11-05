
<a name="getPositions">

#### <span style="color:green">getPositions</span>

Computes the rebalanced and delta numbers of shares for each portfolio
component.

*Call:*

```
getPositions(self, mu, rtype=None, nshares=None, cash=0, ww=None)
```

*Input:*

* `mu` : Rate of reference. Its meaning depends on the optimization method.
For `rtype` set to:
    - `'Risk'` : `mu` is the targeted portfolio expected rate of returns,
    - `'Sharpe'` and `'Sharpe2'` : `mu` is the risk-free rate,
    - `'MinRisk'` and `'InvNRisk'`: `mu` is ignored,
    - `'RiskAverse'` : `mu` is the Lambda aversion coefficient.
* `rtype`: Optimization type. If it is not `None` it will overwrite the value
set by the constructor. The default is `None`.
* `nshares` : Initial number of shares for each portfolio component. The total
value of these shares is the value of the invested capital.
A missing component entry
will be considered `0`. A `None` value assumes that all components entries
are `0`. The name of the components must be present in the `mrkdata`.
The default is `None`.
* `cash` : Additional cash to be added to the capital. A negative entry
assumes a reduction in the total capital  available for rebalance.
The default is `0`.
* `ww` : External portfolio weights (`pd.Series`). If it not set to `None`
these weights will overwrite the calibrated weights.  The default is `None`.

*Returns:* `pd.DataFrame` containing the rolling information.

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
business. Additional cash slippage may occur due to share price differential
between the previous day closing and  execution time.

[TOP](#TOP)

---
