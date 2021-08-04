
<a name="port_view">

#### <span style="color:green">port_view</span>

Plot the portfolio time series together with optional technical indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, view=True, fancy=False)
```

*Input:*

* `emas` :
List for EMA durations. The default is ``[30, 200]``.
* `bollinger` : Boolean flag.
If set `True` it adds the Bollinger bands. The default is `False`.
* `view` : Boolean flag.
`False` suppresses the plotting to tne terminal. The default is `True`.
* `fancy` : Boolean flag with default value 'False'.
    - `False` : it uses the `matplotlib` capabilities.
    - `True` : it uses `plotly` library for interactive time-series view.

*Returns:* pd.DataFrame containing the time-series included in the plot.

[TOP](#TOP)

---
