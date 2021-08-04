
<a name="port_perf">

#### <span style="color:green">port_perf</span>

Brief description of portfolio and its components performances
in terms of average historical rate of returns and maximum drawdowns.

*Call:*

```
port_perf(componly=False, fancy=False)
```

*Input:*

* `componly` : Boolean flag.
If `True`, only the portfolio components maximum drawdowns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format.
    - `True` : The values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame` containing the table of
performance information. Columns:
* `'RR'` : rate of returns
* `'DD'` : maximum rate of drawdown
* `'Beta'` : abs(RR/DD)
* `'DD_date'` : recorder date of maximum drawdown
* `'DD_start'` : start date of maximum drawdown
* `'DD_end'` : end date of maximum drawdown

[TOP](#TOP)

---
