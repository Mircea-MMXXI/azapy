
<a name="port_monthly_returns"></a>

#### <span style="color:green">port_monthly_returns</span>

Computes optimal portfolio and its components monthly (calendar) rate of
returns.

*Call:*

```
port_monthly_returns(withcomp=False, componly=False, fancy=False)
```

*Inputs:*

* `withcomp` : Boolean flag.
If `True`, adds the portfolio components monthly returns to the
report. The default is `False`.
* `componly` : Boolean flag.
If `True`, only the portfolio components monthly returns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format,
    - `True` : The values are reported in percent rounded
    to 2 decimals and presented is color style.

*Returns:* `pd.DataFrame`

[TOP](#TOP)

---
