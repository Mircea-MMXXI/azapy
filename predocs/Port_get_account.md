
<a name="get_account">

#### <span style="color:green">get_account</span>

Returns additional bookkeeping information regarding rebalancing
(*e.g.* residual cash due rounding number of shares, previous period
dividend cash accumulation, etc.)

*Call:*

```
get_account(fancy=False)
```

*Input:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported rounded.

*Returns:* `pd.DataFrame`

Reports, for each rolling period identified by `'Droll'`:

* for each symbol : the number of shares hold,
* `'cash_invst'` : cash invested at the beginning of period,
* `'cash_roll'` : cash rolled to the next period,
* `'cash_divd'` : cash dividend accumulated in the previous period,

> Note: The capital at the beginning of the period is
cash_invst + cash_roll. It is also equal to the previous period:
value of the shares on the fixing date + cash_roll + cash_divd.
There are 2 sources for the cash_roll. The roundup to integer
number of shares and the shares close price differences between
the fixing (computation) and rolling (execution) dates. It could
be positive or negative. The finance of the cash_roll during
each rolling period is assumed  to be done separately by the
investor.

[TOP](#TOP)

---
