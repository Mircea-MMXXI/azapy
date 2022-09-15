
<a name="port_view"></a>

#### <span style="color:green">port_view</span>

Plots the optimal portfolio time series together with some technical
indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, **randomport)
```

*Inputs:*

* `emas` : `list` of int. List of EMA durations. The default is [30, 200].
* `bollinger` : Boolean flag. If it is set `True` then the Bollinger bands are
added. The default is `False`.    
* `opt` : other parameters
    * `fancy` : Boolean flag.
        - `False` : it uses the matplotlib capabilities.
        - `True` : it uses plotly library for interactive time-series view.

        The default is `False`.
    * `title` : `str`. The plot title. The default is `None`.
    * `xaxis` : `str`. The name of x-axis. The default is `'date'`.
    * `yaxis` : `srt`. The name of y-axis. The default is `None`.
    * `saveto` : `str`.
        The name of the file where to save the plot. The default is `None`.

*Returns:* `pandas.DataFrame` containing the time-series included in the plot.

[TOP](#TOP)

---
