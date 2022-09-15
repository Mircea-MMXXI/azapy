#### <span style="color:green">viewFrontiers</span>

Produces a graphical representation of the portfolio frontiers.

*Call:*
```
viewFrontiers(minrisk=True, efficient=20, inefficient=20,
              sharpe=True, musharpe=0,
              component=True, randomport=20, inverseN=True,
              maxdivers=False, divers_efficient=0, divers_inefficient=0,
              addport=None, fig_type='RR_risk', **opt):
```
*Inputs:*
* `minrisk` : Boolean fag. If it is `True` then the minimum risk portfolio will
be visible. The default is `True`.
* `efficient` : Number of points along the optimal frontier (equally spaced
	 along the x-axis). The default is `20`.
* `inefficient` : Number of points along the inefficient frontier (equally
	 spaced along the x-axis). The default is `20`.
* `sharpe` : Boolean flag. If it is `True` then the maximum Sharpe portfolio will
be visible. The default is `True`.
* `musharpe` : Risk-free rate value used in the evaluation of generalized
Sharpe ratio. The default is `0`.
* `component` : Boolean flag. If it is `True` then single asset portfolios
are evaluated and added to the plot for reference. The default is `True`.
* `randomport` : Number of portfolios with random weights (inefficient) to be
evaluate and added to the plot for reference. The default is `20`.
* `inverseN` : Boolean flag. If it is `True` then the equal weighted portfolio and
the optimal portfolio with the same dispersion (risk) value are evaluated and
added to the plot. The default is `True`.
`maxdivers`: Boolean fag. If it is `True` then the maximum diversified portfolio
will be visible. The default is `True`.
* `divers_efficient`: Number of points along the diversified efficient frontier
(equally spaced along the rate of return axis). The default is 20.
* `divers_inefficient`: Number of points along the diversified inefficient frontier
(equally spaced along the rate of return axis). The default is 20.
* `addport` : `dict` or `pandas.DataFrame`.
The weights of additional portfolio to be added to the plot.
If it is a `dict` then the keys are the labels, and the values are
list of weights in the order of `rrate.columns`. If it is a
`pandas.DataFrame` the index are the labels, and each row is a set
of weights. The columns names should match the symbols names.
The default is `None`.
* `fig_type` : `str`. Graphical representation format.
    * `'RR_risk'` : expected rate of return vs risk,
    * `'Sharpe_RR'` : sharpe vs expected rate of return,
    * `'Divers_RR'` : diversification vs expected rate of return.

    The default is `'RR_risk'`.
* `opt` : Additonal parameters:
    * `'title'` : The default is 'Portfolio frontiers'
    * `'xlabel'` : The default is
        - `'risk'` if `fig_type='RR_risk'`
        - `'rate of returns'` otherwise
    * `'ylabel'` : The default is
        - `'rate of returns'` if `fig_type='RR_risk'`
        - `'sharpe'` if `fig_type='RR_sharpe'`
        - `'diversification'` if `fig_type=RR_divers`
    * `'tangent'` : Boolean flag. If set to `True` the tangent
        (to sharpe point) is added. It has effect only  if
        `fig_type='RR_risk'`. The default is `True`.
    * `saveto` : `str`.
        File name to save the figure. The extension dictates the format:
        png, pdf, svg, etc. For more details see the `mathplotlib`
        documentation for `savefig`. The default is `None`.
    * `data` : `dict`.
        Numerical data to construct the plot. If it is not `None` it
        will take precedence and no other numerical evaluations will be
        performed. It is meant to produce different plot representations
        without reevaluations. The default is `None`.

*Returns:* Dictionary containing numerical data used to make the plots.
It can be passed back as `data` argument to reconstruct the plots without
reevaluations.
