
<a name="port_view_all">

#### <span style="color:green">port_view_all</span>

Plot the portfolio and its component time-series in a relative bases.

*Call:*

```
port_view_all(sdate=None, edate=None, componly=False, fancy=False, saveto=None)
```

*Input:*

* `sdate` : `datetime`;
Start date of plotted time-series. If it is set to `None`
then the `sdate` is set to the earliest date in the time-series.
The default is `None`.
* `edate` : `datetime`;
End date of plotted time-series. If it set to `None` then the `edate`
is set to the most recent date of the time-series.
The default is `None`.
* `componly` : Boolean flag with default value `True`.
    - `True` : only the portfolio components time-series are plotted.
    - `False` : the portfolio and its components times-series are plotted.
* `fancy` : Boolean flag with default value `False`.
    - `False` : it uses the `matplotlib` capabilities.
    - `True` : it uses `plotly` library for interactive time-series view.
* `saveto` : File name where to save the plot. The default is `None`.

*Returns:* pd.DataFrame containing the time-series included in the plot.

[TOP](#TOP)

---
