#### <span style="color:green">get_weights</span>

Returns the portfolio weights for each rebalancing period.

*Call:*

```
get_weights(fancy=False)
```

*Inputs:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format.
    - `True` : The values are reported in percent rounded
    to 2 decimals.

*Returns:* `pandas.DataFrame`
