
<a name="getWeights"></a>

#### <span style="color:green">getWeights</span>

Computes the optimal portfolio weights.

*Call:*

```
getWeights(mu, rrate=None, rtype=None, d=1)
```

*Inputs:*

* `mu` : Rate of reference. Its meaning depends on the optimization method.
For `rtype` set to:
    - `'Risk'` : `mu` is the targeted portfolio expected rate of return.
    - `'Sharpe'` and `'Sharpe2'` : `mu` is the risk-free rate.
    - `'MinRisk'` and `'InvNRisk'`: `mu` is ignored.
    - `'RiskAverse'` : `mu` is the risk aversion coefficient $\lambda$.
* `rrate` : `pd.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed in the constructor from `mktdata`. The default is `None`.
* `rtype`: Optimization type. If it is not `None`, it will overwrite the
value set by the constructor. The default is `None`.
* `d` : Frontier type. Has effect only if `rtype='Risk'`. A value of `1` will
trigger the evaluation of optimal portfolio along the efficient frontier.
Otherwise it will find the portfolio with the lowest rate of return along the
inefficient portfolio frontier. The default is `1`.

*Returns:* `pd.Series` containing the portfolio weights.

Note: It will set the following class members:
* _risk_
* _primary_risk_comp_
* _secondary_risk_comp_
* _sharpe_
* _RR_

Their meanings are [here](#RiskMembers).

[TOP](#TOP)

---
