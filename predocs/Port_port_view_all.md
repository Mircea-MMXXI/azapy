#### <span style="color:green">port_view_all</span>

Plots in a relative bases the optimal portfolio and its components time-series.
The components time series prices are designated by the value of
`col_ref` argument in the constructor.

*Call:*

```
port_view_all(sdate=None, edate=None, componly=False, fancy=False, saveto=None)
```

*Inputs:*

* `sdate` : date like;
Start date of plotted time-series. If it is `None`,
then `sdate` is set to the earliest date in the time-series.
The default is `None`.
* `edate` : date like;
End date of plotted time-series. If it is `None`, then `edate`
is set to the most recent date of the time-series.
The default is `None`.
* `componly` : Boolean flag with default value `True`.
    - `True` : only the portfolio components time-series are plotted.
    - `False` : the portfolio and its components times-series are plotted.
* `fancy` : Boolean flag with default value `False`.
    - `False` : it uses the `matplotlib` package capabilities.
    - `True` : it uses `plotly` package for interactive time-series view.
* `saveto` : File name where to save the plot. The extension dictates the
format: `png`, `pdf`, `svg`, etc. For more details see the `mathplotlib`
documentation for `savefig`.The default is `None`.

*Returns:* `pandas.DataFrame` containing the time-series included in the plot.
