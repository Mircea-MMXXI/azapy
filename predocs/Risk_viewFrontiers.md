
<a name="viewFrontiers">

#### <span style="color:green">viewFrontiers</span>

Produces a graphical representation of the portfolio frontiers.

*Call:*
```
viewFrontiers(efficient=20, inefficient=20, musharpe=0., component=True,
              randomport=20, inverseN=True, fig_type='RR_risk',
              options=None, saveto=None, data=None)
```
*Input:*
* `efficient` : Number of points along the optimal frontier (equally spaced
	 along the rate of returns). The default is `20`.
* `inefficient` : Number of points along the inefficient frontier (equally
	 spaced along the rate of returns). The default is `20`.
* `musharpe` : Risk-free rate value used in the evaluation of generalized
Sharpe ratio. The default is `0`.
* `component` : Boolean flag. If `True` the portfolios containing a single
component are evaluated and added to the plot for references.
The default is `True`.
* `randomport` : Number of portfolios with random weights (inefficient) to be
evaluate and added to the plot for reference. The default is `20`.
* `inverseN` : Boolean flag. If `True` the equally weighted portfolio and
the optimal portfolio with the same dispersion (risk) value are evaluated and
added to the plot. The default is `True`.
* `fig_type` : Graphical representation format.   If it is set to `'RR_risk'`
the data is plotted in the rate of return vs dispersion representation,
otherwise the Sharpe vs rate of return will be used. The default is
`'RR_risk'`.
* `options` : A dictionary with additional graphical setups. Relevant keys
are:
    - `'title'` : The default is `'Portfolio frontiers'`.
    - `'xlabel'` : The default is `'risk'` if `fig_type='RR_risk'` and
		 `'rate of returns'` otherwise.
    - `'ylabel'` : The default is `'rate of returns'` if `fig_type='RR_risk'`
		 and `'sharpe'` otherwise.
    - `'tangent'` : Boolean flag. If set to `True` the tangent (to sharpe
		 point) is add. It has effect only  if  `fig_type='RR_risk'`.
		 The default is `True`.
* `saveto` : File name where to save the figure. The extension dictates the
format: `png`, `pdf`, `svg`, etc. For more details see the `mathplotlib`
documentation for `savefig`. The default is `None`.
* `data` : Numerical data to construct the plot. If it is not `None` it
will take precedence and no other numerical evaluations will be
performed. It is meant to produce different plot representations
without reevaluations. The default is `None`.

*Returns:* Dictionary containing numerical data used to make the plots.
It can be passed back to reconstruct the plots without reevaluations.

[TOP](#TOP)

---
