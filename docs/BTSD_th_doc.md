
# BTSD optimal portfolios <a name="TOP"></a>

BTSD stands for Below target Standard Deviation. It is similar
to Delta-risk measure from Omega optimal portfolio but defined in
terms of $L_2$ norm rather than $L_1$,*i.e.*,

\begin{equation*}
  {\rm BTSD}_{\alpha} = \left\| \left( \alpha - r \right)^+ \right\|_2
\end{equation*}

where:

* $\|x\| = ( E[|x|^2])^{1/2}$ is the $L_2$ norm,
* $\alpha$ is the BTSD threshold (it may be interpreted as a risk-free rate),

The BTSD measure can be computed in terms of either standard  or
detrended rate of returns (*i.e.* ${\bar r} = r - E[r]$).

> Note: the BTSD Sharpe ratio for $\alpha=\mu_0$ (the risk-free rate
of return accessible to investor) and standard rate of returns is also known
as Sortino ratio, *i.e.* ${\rm Sortino} = (E[r] - \mu_0)/{\rm BTSD}$.

> Note: BTSD optimal portfolio models with detrended
rate of returns are the same as LSSD first order models.

**azapy** implements a generalization of BTSD measure,
namely the **Mixture BTSD (mBTSD)**.

The mixture is defined as a superposition of regular BTSD measures
for different thresholds, *i.e*,

\begin{equation*}
  \rho = \sum_{l=1}^L {\rm BTSD}_{\alpha_l},
\end{equation*}

where:

* $L$ is the size of the mixture,
* $\{\alpha_l\}_{l=1,\cdots,L}$ is a set of distinct BTSD thresholds.

> Note: a possible choice could be $L=3$ and $\alpha=[0.01, 0.0, -0.01]$

The single BTSD measure is a particular case of mBTSD.

> Note: The mBTSD measures (except for $L=1$, $\alpha_1=0$ and detrended rate
of return) are not proper dispersion measure. They violate the
positive homogeneity axiom and in the case of standard rate of return
they also violate the location invariance axiom.
However, the mathematical formalism of risk-based
optimal portfolio constructions can be applied.

The following portfolio optimization strategies are available:
* Minimization of dispersion for a give expected rate of return,
* Maximization of Sharpe ratio,
* Minimization of the inverse of Sharpe ratio,
* Minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* Maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **BTSDAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_BTSD** : performs portfolio back testing, out-of-sample analysis.

## BTSDAnalyzer class

Computes the portfolio weights and performs in-sample portfolio analysis.

**Methods:**

* [<span style="color:green">getWeights</span>](#getWeights)
* [<span style="color:green">getRsik</span>](#getRisk)
* [<span style="color:green">getPositions</span>](#getPositions)
* [<span style="color:green">viewFrontiers</span>](#viewFrontiers)
* [<span style="color:green">set_mktdata</span>](#set_mktdata)
* [<span style="color:green">set_rrdata</span>](#set_rrdate)
* [<span style="color:green">set_rtype</span>](#set_rtype)
* [<span style="color:green">set_random_seed</span>](#set_random_seed)

Note the following 2 important methods: <a name="RiskMembers"></a>
* **getWeights** : Computes the optimal portfolio weights.
During its computations the following class members are also set:
  * _risk_ : the value of optimal portfolio Delta-risk,
  * _primery_risk_comp_ : redundant (single value list containing the
    optimal portfolio Delta-risk),
  * _secondary_risk_comp_ : redundant (single value list containing the
    optimal portfolio Delta-risk),
  * _sharpe_ : Sharpe ration if `rtype` is set to `'Shapre'` or `'Sharpe2'`
  otherwise `None`. This is the Omega ratio.
  * _RR_ : optimal portfolio expected rate of return.

  * **getPositions** : Provides practical information regarding the portfolio
  rebalancing delta positions and costs.  

### Constructor

```
BTSDAnalyzer(alpha=[0.], coef=None, mktdata=None, colname='adjusted', freq='Q',
             hlength=3.25, calendar=None, rtype='Sharpe', detrended=False,
             method='ecos')
```

where:

* `alpha` : List of distinct BTSD thresholds. The default is `[0.]`.
* `coef` : List of mixture coefficients. Must have the same size as
`alpha`. A `None` value assumes an equal weighted risk mixture.
The default is `None`.
* `mktdata` : `pandas.DataFrame` containing the market data in the format
returned by the function `azapy.readMkT`. The default is `None`.
`mktdata` could be loaded latter.
* `colname` : Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon (portfolio rebalancing period).
It could be `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` years.
* `calendar` :  `np.busdaycalendar` business days calendar. If it is `None`
then the calendar will be set to NYSE business calendar.
The default is `None`.
* `rtype` : Optimization type. The default is `'Sharpe'`. Possible values are:
    - `'Risk'` : minimization of dispersion (risk) measure for a fixed values
    of portfolio expected rate of return,
    - `'Sharpe'` : maximization of generalized Sharpe ratio,
    - `'Sharpe2'` : minimization of inverse generalized Sharpe ratio,
    - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value,
    - `'InvNrisk'` : optimal portfolio with the same dispersion (risk) value
		as equal weighted portfolio,
    - `'RiskAverse'` : optimal portfolio for a fixed risk aversion coefficient.
* `detrended` : Boolean flag.
In the BTSD expression use:
    - `True` : detrended rate of return, *i.e.* ${\bar r} = r - E[r]$,
    - `False` : standard rate of return.
The default is `False`.
* `method` : Designates the SOCP numerical method.
It could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.

> Note:
>	* `'ecos'` : is SOCP implementation of **ecos** *(Embedded Cone Solver)*
package.
> * `'cvxopt'` : is the SOCP implantation from **cvxopt** package.
>
> In our cases `'ecos'` is the fastest.

[TOP](#TOP)

### Methods:

<a name="getWeights"></a>

#### <span style="color:green">getWeights</span>

Computes the optimal portfolio weights.

*Call:*

```
getWeights(mu, rrate=None, rtype=None, d=1)
```

*Inputs:*

* `mu` : Rate of reference. Its meaning depends on the optimization method.
For `rtype` set to:
    - `'Risk'` : `mu` is the targeted portfolio expected rate of return.
    - `'Sharpe'` and `'Sharpe2'` : `mu` is the risk-free rate.
    - `'MinRisk'` and `'InvNRisk'`: `mu` is ignored.
    - `'RiskAverse'` : `mu` is the risk aversion coefficient $\lambda$.
* `rrate` : `pd.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed in the constructor from `mktdata`. The default is `None`.
* `rtype`: Optimization type. If it is not `None`, it will overwrite the
value set by the constructor. The default is `None`.
* `d` : Frontier type. Has effect only if `rtype='Risk'`. A value of `1` will
trigger the evaluation of optimal portfolio along the efficient frontier.
Otherwise it will find the portfolio with the lowest rate of return along the
inefficient portfolio frontier. The default is `1`.

*Returns:* `pd.Series` containing the portfolio weights.

Note: It will set the following class members:
* _risk_
* _primary_risk_comp_
* _secondary_risk_comp_
* _sharpe_
* _RR_

Their meanings are [here](#RiskMembers).

[TOP](#TOP)

---

<a name="getRisk"></a>

#### <span style="color:green">getRisk</span>

Computes the risk of a portfolio defined by a set of weights.

*Call:*
```
getRisk(ww, rrate=None)
```

*Inputs:*

* `ww` : List like of portfolio weights. Its length must be equal to the
number of symbols in `rrate` (`mktdata`). All weights must by $\ge 0$ and
their sum equal to $1$. If it
is a `list` or a `np.array` then the weights are assumed to be in the order
of `rrate.columns`. If it is a `pd.Series` the index should be compatible
with the `rrate.columns` or `mktdata` symbols (not necessary in the same
order).
* `'rrate'` : `pd.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed by the constructor from `mktdata`. The default is `None`.

*Returns:* The value of the risk measure.

Note: It will set the following class members:
* _risk_
* _primary_risk_comp_
* _secondary_risk_comp_
* _sharpe_
* _RR_

Their meanings are [here](#RiskMembers).

[TOP](#TOP)

---

<a name="getPositions"></a>

#### <span style="color:green">getPositions</span>

Computes the rebalanced and delta numbers of shares for each portfolio
component.

*Call:*

```
getPositions(mu, rtype=None, nshares=None, cash=0, ww=None)
```

*Inputs:*

* `mu` : Rate of reference. Its meaning depends on the optimization method.
For `rtype` set to:
    - `'Risk'` : `mu` is the targeted portfolio expected rate of returns,
    - `'Sharpe'` and `'Sharpe2'` : `mu` is the risk-free rate,
    - `'MinRisk'` and `'InvNRisk'`: `mu` is ignored,
    - `'RiskAverse'` : `mu` is the risk aversion coefficient $\lambda$.
* `rtype`: Optimization type. If it is not `None`, it will overwrite the value
set by the constructor. The default is `None`.
* `nshares` : Initial number of shares for each portfolio component. The total
value of these shares is the value of the invested capital.
A missing component entry
will be considered `0`. A `None` value assumes that all components entries
are `0`. The name of the components must be present in the `mrkdata`.
The default is `None`.
* `cash` : Additional cash to be added to the capital. A negative entry
assumes a reduction in the total capital  available for rebalance.
The default is `0`.
* `ww` : `pd.Series` external portfolio weights. If it is not `None`
these weights will overwrite the calibrated weights. The default is `None`.

*Returns:* `pd.DataFrame` containing the rolling information.

Columns:

* `'old_nsh'` : initial number of shares per portfolio component as well as
additional cash position. These are present in the input.
* `'new_nsh'` : the new number of shares per component plus the  residual
cash (due to the rounding to an integer number of shares). A negative entry
means that the investor needs to add more cash in order to cover for the
number of share  roundups. It has a small value.
* `'diff_nsh'` : delta number of shares - the number of shares that needs to be
both/sold in order to rebalance the portfolio positions.
* `'weights'` : portfolio weights used for rebalance. The `'cash'` entry
is the new portfolio value (invested capital).
* `'prices'` : share prices used for rebalance evaluations.

>Note: Since the prices are closing prices, the rebalance can be executed next
business day. Additional cash slippage may occur due to share price differential
between the previous day closing and  execution time.

[TOP](#TOP)

---

<a name="viewFrontiers"></a>

#### <span style="color:green">viewFrontiers</span>

Produces a graphical representation of the portfolio frontiers.

*Call:*
```
viewFrontiers(efficient=20, inefficient=20, musharpe=0.,
              component=True, randomport=20, inverseN=True,
              fig_type='RR_risk', options=None, saveto=None,
              data=None)
```
*Inputs:*
* `efficient` : Number of points along the optimal frontier (equally spaced
	 along the x-axis). The default is `20`.
* `inefficient` : Number of points along the inefficient frontier (equally
	 spaced along the x-axis). The default is `20`.
* `musharpe` : Risk-free rate value used in the evaluation of generalized
Sharpe ratio. The default is `0`.
* `component` : Boolean flag. If `True` the portfolios containing a single
component are evaluated and added to the plot for reference.
The default is `True`.
* `randomport` : Number of portfolios with random weights (inefficient) to be
evaluate and added to the plot for reference. The default is `20`.
* `inverseN` : Boolean flag. If `True` the equal weighted portfolio and
the optimal portfolio with the same dispersion (risk) value are evaluated and
added to the plot. The default is `True`.
* `fig_type` : Graphical representation format. If it is set to `'RR_risk'`
the data is plotted in the risk vs rate of return representation,
otherwise the rate of return vs Sharpe will be used. The default is
`'RR_risk'`.
* `options` : A dictionary with additional graphical setups. Relevant keys
are:
    - `'title'` : The default is `'Portfolio frontiers'`.
    - `'xlabel'` : The default is `'risk'` if `fig_type='RR_risk'` and
		 `'rate of returns'` otherwise.
    - `'ylabel'` : The default is `'rate of returns'` if `fig_type='RR_risk'`
		 and `'sharpe'` otherwise.
    - `'tangent'` : Boolean flag. If set to `True` the tangent (to sharpe
		 point) is added. It has effect only  if  `fig_type='RR_risk'`.
		 The default is `True`.
* `saveto` : File name where to save the figure. The extension dictates the
format: `png`, `pdf`, `svg`, etc. For more details see the `mathplotlib`
documentation for `savefig`. The default is `None`.
* `data` : Precomputed numerical data used to construct the plot.
If it is not `None` it
will take precedence and no other numerical evaluations will be
performed. Its main use is to produce different plot representations
without reevaluations. The default is `None`.

*Returns:* Dictionary containing numerical data used to make the plots.
It can be passed back as `data` argument to reconstruct the plots without
reevaluations.

[TOP](#TOP)

---

<a name="set_mktdata"></a>

#### <span style="color:green">set_mktdata</span>

Sets historical market data. It will overwrite the choices made in the
constructor.

*Call:*

```
set_mktdata(mktdata, colname='adjusted', freq='Q', hlength=3.25, calendar=None)
```

*Inputs:*

* `mktdata` : `pd.DataFrame`
Historic daily market data for portfolio components in the format
returned by `azapy.mktData` function.
* `colname` :
Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` :
Rate of returns horizon. It could be
`'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` :
History length in number of years used for calibration. A
fractional number will be rounded to an integer number of months.
The default is `3.25` years.
* `calendar` : `np.busdaycalendar`, optional
Business days calendar. If it is `None`, then the calendar will be set
to NYSE business calendar. The default is `None`.


*Returns:* `None`

[TOP](#TOP)

---

<a name="set_rrate"></a>

#### <span style="color:green">set_rrate</span>

Sets portfolio components historical rates of returns.
It will overwrite the value computed by the constructor from `mktdata`.

*Call:*

```
set_rrate(rrate)
```

*Inputs:*

* `rrate` : `pd.DataFrame`,
portfolio components historical rates of returns, where the
columns are `'date'`, `symbol1`, `symbol2`, etc.


*Returns:* `None`

[TOP](#TOP)

---

<a name="set_rtype"></a>

#### <span style="color:green">set_rtype</span>

Sets the optimization type. It will overwrite the value set in the
constructor.

*Call:*

```
set_rtype(rtype)
```

*Inputs:*

* `rtype` : Optimization type.

*Returns:* `None`

[TOP](#TOP)

---

<a name="set_random_seed"></a>

#### <span style="color:green">set_random_seed</span>

Sets the seed for Dirichlet random generator used in `viewFrontiers`.

*Call:*

```
set_random_seed(seed=42)
```

*Inputs:*

* `seed` : The random generator seed - in case you want to set it to a weird
value other than 42 :). The default is `42`.

*Returns:* `None`

[TOP](#TOP)

---
<a name="OmegaAnalyzer_class_example"></a>

### [Examples](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/analyzers/BTSDAnalyzer_examples.py)

```
import numpy as np
import pandas as pd
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Set the BTSD mixture parameter
alpha = [0.01, 0, -0.01]
coef = [1, 1, 2]

#=============================================================================
# Compute Sharpe optimal portfolio
# build the analyzer object
cr1 = az.BTSDAnalyzer(alpha, coef, mktdata)
# computes Sharpe weights for 0 risk-free rate
ww1 = cr1.getWeights(mu=0.)
# print portfolio characteristics
# primary risk = [Delta-risk] (redundant)
# secondary risk = [Delta-risk] (redundant)
# risk = Delta-risk
# Share = BTSD-Sharpe ratio
RR = cr1.RR
risk = cr1.risk
prim = cr1.primary_risk_comp.copy()
seco = cr1.secondary_risk_comp.copy()
sharpe = cr1.sharpe
print("\nSharpe optimal portfolio\n")
print(f"status {cr1.status}")
print(f"coef {ww1}")
print(f"Secondary risk {seco}")
print(f"Primary risk {prim}")
print(f"Sharpe {sharpe}")
print(f"RR {RR}")
print(f"risk {risk}")

# Test risk by computing the risk of a portfolio with weights ww1
test_risk = cr1.getRisk(ww1)
test_risk_res = pd.DataFrame({'risk': [risk], 'test_risk': [test_risk],
                              'diff': [risk-test_risk]})
print(f"Test for the risk computation\n {test_risk_res}")

# Test the Sharpe weights by estimating an optimal portfolio with
# the same expected rate of returns.
test_ww1 = cr1.getWeights(mu=RR, rtype='Risk')
ww_comp = pd.DataFrame({"ww1": ww1, "test_ww1": test_ww1,
                        'diff': ww1-test_ww1})
print(f"Test for weights computation\n {ww_comp}")

#=============================================================================
# Frontiers evaluations
print("\nFrontiers evaluations\n")
opt = {'title': "BTSD Port", 'tangent': True}
print("\n rate of returns vs risk representation")
rft = cr1.viewFrontiers(musharpe=0, randomport=100, options=opt)
print("\n Sharpe vs rate of returns representation")
rft2 = cr1.viewFrontiers(data=rft, fig_type='Sharpe_RR')

#=============================================================================
# Sharpe vs. Sharpe2
# first Sharpe (default rtype)
cr1 = az.BTSDAnalyzer(alpha, coef, mktdata)
ww1 = cr1.getWeights(mu=0.)
RR1 = cr1.RR
risk1 = cr1.risk
prim1 = cr1.primary_risk_comp.copy()
seco1 = cr1.secondary_risk_comp.copy()
sharpe1 = cr1.sharpe
# second Sharpe2
cr2 = az.BTSDAnalyzer(alpha, coef, mktdata)
ww2 = cr2.getWeights(mu=0., rtype="Sharpe2")
RR2 = cr2.RR
risk2 = cr2.risk
prim2 = cr2.primary_risk_comp.copy()
seco2 = cr2.secondary_risk_comp.copy()
sharpe2 = cr2.sharpe
# print comparison - must be very close
print("\nSharpe vs. Sharpe2\n")
print(f"status {cr2.status} = {cr1.status}")
ww_comp = pd.DataFrame({"ww2": ww2, "ww1": ww1, "diff": ww2-ww1})
print(f"coef\n {ww_comp}")
seco_comp = pd.DataFrame({"seco2": seco2, "seco1": seco1, "diff": seco2-seco1})
print(f"Secondary risk\n {seco_comp}")
prim_comp = pd.DataFrame({"prim2": prim2, "prim1": prim1,
                          "diff": prim2-prim1})
print(f"Primary risk\n {prim_comp}")
RR_comp = pd.DataFrame({'RR2': [RR2], 'RR1': [RR1], 'diff': [RR2 - RR1]})
print(f"RR comp\n {RR_comp}")
risk_comp = pd.DataFrame({'risk2': [risk2], 'risk1': [risk1],
                          'diff': [risk2-risk1]})
print(f"risk comp\n {risk_comp}")
sharpe_comp = pd.DataFrame({'sharpe2': [sharpe2], 'sharpe1': [sharpe1],
                            'diff': [sharpe2-sharpe1]})
print(f"Sharpe comp\n {sharpe_comp}")

# # Speed of Sharpe vs Sharpe2 - may take some time
# # please uncomment the lines below
# %timeit cr2.getWeights(mu=0., rtype='Sharpe')
# %timeit cr2.getWeights(mu=0., rtype='Sharpe2')

#=============================================================================

# Compute InvNrisk optimal portfolio
cr1 = az.BTSDAnalyzer(alpha, coef, mktdata)
# compute the weights of InvNrisk
ww1 = cr1.getWeights(mu=0., rtype="InvNrisk")
RR1 = cr1.RR

# Test - compute the optimal portfolio for RR1 targeted rate of return
ww2 = cr1.getWeights(mu=RR1, rtype="Risk")
# print comparison results - must be very close
print("\nInvNrisk\n")
ww_comp = pd.DataFrame({"InvNrisk": ww1, "Optimal": ww2, 'diff': ww1-ww2})
print(f"weights comp\n {ww_comp}")

# Test - compute the risk of equal weighted portfolio
ww = np.ones(len(symb))
ww = ww / np.sum(ww)
risk = cr1.getRisk(ww)
# print comparison results - must be identical
risk_comp = pd.DataFrame({'1/N': [risk], 'InvNrisk': [cr1.risk],
                          'diff': [risk - cr1.risk]})
print(f"risk comp\n {risk_comp}")

#=============================================================================
# Compute MinRisk optimal portfolio
cr1 = az.BTSDAnalyzer(alpha, coef, mktdata)
# compute the MinRisk portfolio
ww1 = cr1.getWeights(mu=0., rtype="MinRisk")

# Test - using rtype='Risk' for expected rate of return 0
# should default to 'MinRisk' optimal portfolio
ww2 = cr1.getWeights(mu=0., rtype="Risk")
# print comparison - should be identical
print("\nMinRisk\n")
ww_comp = pd.DataFrame({"MinRisk": ww1, "Test": ww2, 'diff': ww1-ww2})
print(f"weights comp\n {ww_comp}")

#=============================================================================
# Compute RiskAverse optimal portfolio
# first compute the Sharpe portfolio
cr1 = az.BTSDAnalyzer(alpha, coef, mktdata)
ww1 = cr1.getWeights(mu=0.)
sharpe = cr1.sharpe
risk = cr1.risk

# compute RiskAverse portfolio for Lambda=sharpe
Lambda = sharpe
cr2 = az.BTSDAnalyzer(alpha, coef, mktdata)
ww2 = cr2.getWeights(mu=Lambda, rtype='RiskAverse')

# comparison - they should be very close
print("\nRiskAverse\n")
risk_comp = pd.DataFrame({'risk': [cr2.risk], 'test': [cr2.RR / Lambda],
                          'Sharpe risk': [risk]})
print(f"risk comp\n {risk_comp}")
ww_comp = pd.DataFrame({'ww1': ww1, 'ww2': ww2, 'diff': ww1-ww2})
print(f"weigths:\n {ww_comp}")

#=============================================================================
# # speed comparisons for different SOCP methods
# # may take some time to complete
# # please uncomment the lines below
# import time
# methods = ['ecos', 'cvxopt']
# xta = {}
# for method in methods:
#     crrx = az.BTSDAnalyzer(alpha, coef, mktdata, method=method)
#     toc = time.perf_counter()
#     wwx = crrx.getWeights(mu=0.)
#     tic = time.perf_counter() - toc
#     print(f"method: {method} time: {tic}")
#     xta[method] = pd.Series([tic], index=["Time"]).append(wwx)

# res = pd.DataFrame(xta)
# print(res.round(4))

#=============================================================================
# Example of rebalancing positions
cr1 = az.BTSDAnalyzer(alpha, coef, mktdata)

# existing positions and cash
ns = pd.Series(100, index=symb)
cash = 0.

# new positions and rolling info
pos = cr1.getPositions(mu=0., rtype='Sharpe', nshares=ns, cash=0.)
print(f" New position report\n {pos}")
```

[TOP](#TOP)

---

## Port_BTSD class

Out-of-Sample (back testing) simulation of mBTSD optimal portfolio periodically
rebalanced.


**Methods:**

* [<span style="color:green">set_model</span>](#set_model)
* [<span style="color:green">port_view</span>](#port_view)
* [<span style="color:green">port_view_all</span>](#port_view_all)
* [<span style="color:green">port_drawdown</span>](#port_drawdown)
* [<span style="color:green">port_perf</span>](#port_perf)
* [<span style="color:green">port_annual_returns</span>](#port_annual_returns)
* [<span style="color:green">port_monthly_returns</span>](#port_monthly_returns)
* [<span style="color:green">port_period_returns</span>](#port_period_returns)
* [<span style="color:green">get_nshares</span>](#get_nshares)
* [<span style="color:green">get_weights</span>](#get_weights)
* [<span style="color:green">get_account</span>](#get_account)
* [<span style="color:green">get_mktdata</span>](#get_mktdata)


The most important method is **set_model**. It must be called before any
other method.

### Constructor

```
Port_BTSD(mktdata, symb=None, sdate=None, edate=None, col_price='close',
          col_divd='divd', col_ref='adjusted', col_calib='adjusted',
          pname='Port', pcolname=None, capital=100000, schedule=None,
          freq='Q', noffset=-3, fixoffset=-1, calendar=None)
```


where:

* `mktdata` : `pd.DataFrame`;
Market data in the format `"symbol"`, `"date"`, `"open"`, `"high"`,
`"low"`, `"close"`, `"volume"`, `"adjusted"`, `"divd"`, `"split"`
(*e.g.* as returned by `azapy.readMkT`).
* `symb` :
List of symbols of portfolio components. All symbols
should be present in `mktdata`. If it is `None`, then `symb` will default
to the full set of symbols present in `mktdata`. The default
is `None`.
* `sdate` : date like;
Start date for historical simulation. If it is `None`, then `sdate` will
default to the earliest date in `mktdata`. The default is `None`.
* `edate` : date like;
End date for historical simulation. Must be
greater than  `sdate`. If it is `None`, then `edate` will default
to the latest date in `mktdata`. The default is `None`.
* `col_price` : `string`;
Column name in the `mktdata` that will be considered
for portfolio aggregation. The default is `'close'`.
* `col_divd` : `string`;
Column name in the `mktdata` that holds the dividend
information. The default is `'dvid'`.
* `col_ref` : `string`;
Column name in the `mktdata` that will be used as a price
reference for portfolio components (used for various comparisons and graphs).
The default is `'adjusted'`.
* `col_calib` : `string`;
Column name in the `mktdata` used for historical weights calibrations.
The default is `'adjusted'`.
* `pname` : `string`;
The name of the portfolio. The default is `'Port'`.
* `pcolname` : `string`;
Name of the portfolio price column. If it is `None`, than
`pcolname=pname`. The default is `None`.
* `capital` : `float`;
Initial portfolio Capital in dollars. The default is `100000`.
* `schedule` : `pd.DataFrame`;
Rebalancing schedule, with columns for `'Droll'` rolling date and
`'Dfix'` fixing date. If it is `None` than the schedule will be set
using the `freq`, `nsoffset`, `fixoffset`, `hlength` and `calendar`
information. The default is `None`.
* `freq` : `string`;
Rebalancing frequency. It can be `'Q'` for quarterly or `'M'` for
monthly rebalancing. It is relevant only if schedule
is `None`. The default is `'Q'`.
* `noffset` : `int`;
Rebalancing date `'Droll'` number of offset business days
relative to the end of the period (quart or month). A positive
value add business days beyond the calendar end of the period while
a negative value subtract business days. It is relevant only if
`schedule` is `None`. The default is `-3`.
* `fixoffset` : `int`;
fixing date `'Dfix'` number of offset business days relative to
the rebalancing date `'Droll'`. It cane be `0` or negative. It is
relevant only if `schedule` is `None`. The default is `-1`.
* `calendar` : `np.busdaycalendar`;
Business calendar. If it is `None`, then it will be set to NYSE
business calendar. The default is `None`.

[TOP](#TOP)

### Methods:

<a name="set_model">

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
set_model(mu, alpha=[0.], coef=None, rtype='Sharpe', detrended=False,
          hlength=3.25, method='ecos'):
```

*Inputs:*

* `mu` :
Reference rate. Its meaning depends of the value of `rtype`. For
`rtype` equal to:
    - `'Risk'` : `mu` is the targeted expected rate of returns,
    - `'Sharpe'` and `'Sharpe2'`: `mu` is the risk-free rate,
    - `'MinRisk'` and `'InvNrisk'` : `mu` is ignored,
    - `'RiskAverse'` : `mu` is the risk aversion coefficient $\lambda$.
* `alpha` : List of distinct BTSD thresholds. The default is `[0.]`.
* `coef` : List of mixture coefficients. Must have the same size as
`alpha`. A `None` value assumes an equal weighted risk mixture.
The default is `None`.
* `rtype` :
Optimization type. The default is `'Sharpe'`. Possible values are:
    - `'Risk'` : minimization of dispersion (risk) measure for a fixed values
    of portfolio expected rate of return,
    - `'Sharpe'` : maximization of generalized Sharpe ratio,
    - `'Sharpe2'` : minimization of inverse generalized Sharpe ratio,
    - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value,
    - `'InvNrisk'` : optimal portfolio with the same dispersion (risk) value
		as equal weighted portfolio,
    - `'RiskAverse'` : optimal portfolio for a fixed risk aversion coefficient.
* `detrended` : Boolean flag.
In the BTSD expression use:
    - `True` : detrended rate of return, *i.e.* ${\bar r} = r - E[r]$,
    - `False` : standard rate of return.
The default is `False`.
* `hlength` :
The length in years of historical calibration period relative
to `'Dfix'`. A fractional number will be rounded to an integer number
of months. The default is `3.25` years.
* `method` :
Designates the SOCP method.
It could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---

<a name="port_view"></a>

#### <span style="color:green">port_view</span>

Plots the optimal portfolio time series together with some technical
indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, fancy=False, saveto=None)
```

*Inputs:*

* `emas` :
List for EMA durations. The default is ``[30, 200]``.
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

*Returns:* `pd.DataFrame` containing the time-series included in the plot.

[TOP](#TOP)

---

<a name="port_view_all"></a>

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

*Returns:* `pd.DataFrame` containing the time-series included in the plot.

[TOP](#TOP)

---

<a name="port_drawdown"></a>

#### <span style="color:green">port_drawdown</span>

Computes the portfolio drawdowns.

*Call:*

```
port_drawdown(top=5, fancy=False)
```

*Inputs:*

* `top` :
The number of largest drawdown that will be reported.
The default is `5`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame` containing the table of
drawdown events. Columns:
* `'DD'` : drawdown rate,
* `'Date'` : recorded date of the drawdown,
* `'Star'` : start date of the drawdown,
* `'End'` : end date of the drawdown. A `NaN` value indicates that the
drawdown event is in progress and the values of `'DD'` and `'Date'` are
provisional only.

[TOP](#TOP)

---

<a name="port_perf"></a>

#### <span style="color:green">port_perf</span>

Brief description of optimal portfolio and its components performances
in terms of average historical rate of returns and maximum drawdowns.

*Call:*

```
port_perf(componly=False, fancy=False)
```

*Inputs:*

* `componly` : Boolean flag.
If `True`, only the portfolio components information is reported.
The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame` containing the table of
performance information. Columns:
* `'RR'` : annual average rate of returns,
* `'DD'` : maximum rate of drawdown during the simulation period,
* `'Beta'` : `abs(RR/DD)`,
* `'DD_date'` : recorded date of maximum drawdown,
* `'DD_start'` : start date of maximum drawdown,
* `'DD_end'` : end date of maximum drawdown.

[TOP](#TOP)

---

<a name="port_annual_returns"></a>

#### <span style="color:green">port_annual_returns</span>

Computes optimal portfolio and its components annual (calendar) rates of returns.
The components time series prices used in the estimations are designated by
the value of `col_ref` argument in the constructor.

*Call:*

```
port_annual_returns(withcomp=False, componly=False, fancy=False)
```

*Inputs:*

* `withcomp` : Boolean flag.
If `True`, adds the portfolio components annual returns to the
report. The default is `False`.
* `componly` : Boolean flag.
If `True`, only the portfolio components annual returns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals and presented is color style.

*Returns:* `pd.DataFrame`

[TOP](#TOP)

---

<a name="port_monthly_returns"></a>

#### <span style="color:green">port_monthly_returns</span>

Computes optimal portfolio and its components monthly (calendar) rate of
returns.

*Call:*

```
port_monthly_returns(withcomp=False, componly=False, fancy=False)
```

*Inputs:*

* `withcomp` : Boolean flag.
If `True`, adds the portfolio components monthly returns to the
report. The default is `False`.
* `componly` : Boolean flag.
If `True`, only the portfolio components monthly returns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals and presented is color style.

*Returns:* `pd.DataFrame`

[TOP](#TOP)

---

<a name="port_period_returns"></a>

#### <span style="color:green">port_period_returns</span>

Computes the rolling periods rate of returns.

*Call:*

```
port_period_returns(fancy=False)
```

*Inputs:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame`

Each rolling period is indicated by its start date, `Droll`.
For reference, the values of `Dfix` and components weights are
included in the report.

[TOP](#TOP)

---

<a name="get_nshares"></a>

#### <span style="color:green">get_nshares</span>

Returns the number of shares hold after each rolling date.

*Call:*

```
get_nshares()
```

*Inputs:* None


*Returns:* `pd.DataFrame`

Each rolling period is indicated by its start date, `Droll`.


[TOP](#TOP)

---

<a name="get_account"></a>

#### <span style="color:green">get_account</span>

Returns additional bookkeeping information regarding rebalancing
(*e.g.* residual cash due the number of shares roundup to an integer,
previous period dividend cash accumulation, etc.)

*Call:*

```
get_account(fancy=False)
```

*Inputs:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported rounded.

*Returns:* `pd.DataFrame`

Accounting report; each rolling period is identified by `'Droll'`. Columns:

* for each symbol : number of shares hold,
* `'cash_invst'` : cash invested at the beginning of the period,
* `'cash_roll'` : cash rolled to the next period,
* `'cash_divd'` : cash dividend accumulated in the previous period.

> Note: The capital at the beginning of the rolling period is
`'cash_invst'` + `'cash_roll'`. It is also equal to the previous period
value of the shares on the fixing date + `'cash_roll'` + `'cash_divd'`.
There are 2 sources for `'cash_roll'`. The roundup to an integer
number of shares and the shares price differential between
the fixing (computation) and rolling (execution) dates. In general it
has a small positive or negative value.
The finance of the `'cash_roll'` (if it has a negative value) is assumed
to be done separately by the investor.

[TOP](#TOP)

---

<a name="get_mktdata"></a>

#### <span style="color:green">get_mktdata</span>

Returns the actual market data used for portfolio evaluations.

*Call:*

```
get_mktdata()
```

*Inputs:* None


*Returns:* `pd.DataFrame`

[TOP](#TOP)

---

### [Examples](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_BTSD_examples.py)

```
import time
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute BTSD-Sharpe optimal portfolio
alpha = [0.01, 0., -0.01]
coef = [1, 2, 3]

p4 = az.Port_BTSD(mktdata, pname='BTSDPort')

tic = time.perf_counter()
port4 = p4.set_model(mu=0., alpha=alpha, coef=coef)   
toc = time.perf_counter()
print(f"time Sharpe: {toc-tic}")

ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)

# Use rtype='Sharpe2' - should be the same results
tic = time.perf_counter()
port4_2 = p4.set_model(mu=0., alpha=alpha, coef=coef, rtype='Sharpe2')   
toc = time.perf_counter()
print(f"time Sharpe2: {toc-tic}")

# compare - should be identical
port4.columns = ['Sharpe']
port4_2.columns = ['Sharpe2']
pp = az.Port_Simple([port4, port4_2])
_ = pp.set_model()
_ = pp.port_view_all(componly=(True))

#=============================================================================
# Compute BTSD optimal portfolio
port4 = p4.set_model(mu=0.1, alpha=alpha, coef=coef, rtype="Risk")   
ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)

#=============================================================================
# Compute minimum BTSD optimal portfolio
port4 = p4.set_model(mu=0.1, alpha=alpha, coef=coef, rtype="MinRisk")   
ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)

#=============================================================================
# Compute optimal portfolio with BTSD of equal weighted portfolio
port4 = p4.set_model(mu=0.1, alpha=alpha, coef=coef, rtype="InvNrisk")   
ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)

#=============================================================================
# Compute optimal portfolio for fixed risk aversion
port4 = p4.set_model(mu=0.5, alpha=alpha, coef=coef, rtype="RiskAverse")  
ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)  

#=============================================================================
# # speed comparisons for different LP methods
# # may take some time to complete
# # please uncomment the lines below
# methods = ['ecos', 'highs-ds', 'highs-ipm', 'highs', 'glpk', 'cvxopt',  
#            'interior-point' ]
# zts = []
# for method in methods:
#     toc = time.perf_counter()
#     zz = p4.set_model(mu=0., alpha=alpha, coef=coef, method=method)  
#     tic = time.perf_counter()
#     print(f"{method} time: {tic-toc}")  
#     zz.columns = [method]
#     zts.append(zz)

# # must be identical   
# pp = az.Port_Simple(zts)
# _ = pp.set_model()
# _ = pp.port_view_all(componly=True)
```

[TOP](#TOP)