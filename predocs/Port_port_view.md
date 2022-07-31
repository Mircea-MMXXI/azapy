#### <span style="color:green">port_view</span>

Plots the optimal portfolio time series together with some technical
indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, fancy=False, saveto=None)
```

*Inputs:*

* `emas` :
List for EMA durations. The default is `[30, 200]`.
* `bollinger` : Boolean flag.
`True` adds the Bollinger bands. The default is `False`.
* `view` : Boolean flag.
`False` suppresses the plotting to the terminal. The default is `True`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : it uses the `matplotlib` package capabilities.
    - `True` : it uses `plotly` package for interactive time-series view.
* `saveto` : File name where to save the plot. The extension dictates the
format: `png`, `pdf`, `svg`, etc. For more details see the `mathplotlib`
documentation for `savefig`. The default is `None`.

*Returns:* `pandas.DataFrame` containing the time-series included in the plot.
