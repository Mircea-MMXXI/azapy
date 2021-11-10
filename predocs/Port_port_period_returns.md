
<a name="port_period_returns"></a>

#### <span style="color:green">port_period_returns</span>

Computes the rolling periods rate of returns.

*Call:*

```
port_period_returns(fancy=False)
```

*Inputs:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format.
    - `True` : The values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame`

Each rolling period is indicated by its start date, `Droll`.
For reference, the values of `Dfix` and components weights are
included in the report.

[TOP](#TOP)

---
