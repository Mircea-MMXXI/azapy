#### <span style="color:green">port_view_all</span>

Plots in a relative bases the optimal portfolio and its components time-series.
The components time series prices are designated by the value of
`col_ref` argument in the constructor.

*Call:*

```
port_view_all(sdate=None, edate=None, componly=False, **opt)
```

*Inputs:*

* `sdate` : date-like;
Start date of plotted time-series. If it is `None`,
then `sdate` is set to the earliest date in the time-series.
The default is `None`.
* `edate` : date-like;
End date of plotted time-series. If it is `None`, then `edate`
is set to the most recent date of the time-series.
The default is `None`.
* `componly` : Boolean flag.
    - `True` : only the portfolio components time-series are plotted.
    - `False`: the portfolio and its components times-series are plotted.

    The default is `True`.
* `opt` : Other parameters:
    * `fancy` : Boolean flag.
        - `False` : it uses the pandas plot (matplotlib) capabilities.
        - `True` : it uses plotly library for interactive time-series view.

        The default is `False`.
    * `title` : `str`. The plot title. The default is `None`.
    * `xaxis` : `str`. The name of x-axis. The default is `'date'`.
    * `yaxis` : `srt`. The name of y-axis. The default is `None`.
    * `saveto` : `str`.
        The name of the file where to save the plot. The default is `None`.

*Returns:* `pandas.DataFrame` containing the time-series included in the plot.
