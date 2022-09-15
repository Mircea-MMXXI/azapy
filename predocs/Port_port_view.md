#### <span style="color:green">port_view</span>

Plots the optimal portfolio time series together with some technical
indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, **randomport)
```

*Inputs:*

* `emas` : `list` of int, optional; List of EMA durations. The default is [30, 200].
* `bollinger` : `Boolean`, optional; If it is `True` then the Bollinger bands are
added. The default is `False`.    
* `opt` : other parameters
    * `fancy` : `Boolean`, optional;
        - `False` : it uses the matplotlib capabilities.
        - `True` : it uses plotly library for interactive time-series view.

        The default is `False`.
    * `title` : `str`, optional; The plot title. The default is `'Relative performance'`.
    * `xaxis` : `str`, optional; The name of x-axis. The default is `'date'`.
    * `yaxis` : `str`, optional; The name of y-axis. The default is `None`.
    * `saveto` : `str`, optional;
        The name of the file where to save the plot. The default is `None`.

*Returns:* `pandas.DataFrame` containing the time-series included in the plot.
