#### <span style="color:green">get_account</span>

Returns additional bookkeeping information regarding rebalancing
(*e.g.* residual cash due the number of shares roundup to an integer,
previous period dividend cash accumulation, etc.)

*Call:*

```
get_account(fancy=False)
```

*Inputs:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported rounded.

*Returns:* `pd.DataFrame`

Accounting report; each rolling period is identified by `'Droll'`. Columns:

* for each symbol : number of shares hold,
* `'cash_invst'` : cash invested at the beginning of the period,
* `'cash_roll'` : cash rolled to the next period,
* `'cash_divd'` : cash dividend accumulated in the previous period.

> Note: The capital at the beginning of the rolling period is
`'cash_invst'` + `'cash_roll'`. It is also equal to the previous period
value of the shares on the fixing date + `'cash_roll'` + `'cash_divd'`.
There are 2 sources for `'cash_roll'`. The roundup to an integer
number of shares and the shares price differential between
the fixing (computation) and rolling (execution) dates. In general it
has a small positive or negative value.
The finance of the `'cash_roll'` (if it has a negative value) is assumed
to be done separately by the investor.
