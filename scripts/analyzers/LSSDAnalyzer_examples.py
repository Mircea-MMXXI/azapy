# Examples
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
# Define mLSSD measure parameters coef
coef = np.ones(3)
coef = coef / coef.sum()

#=============================================================================
# Compute Sharpe optimal portfolio
# build the analyzer object
cr1 = az.LSSDAnalyzer(coef, mktdata)
# computes Sharpe weights for 0 risk-free rate
ww1 = cr1.getWeights(mu=0.)
# print portfolio characteristics
# primary risk = set of LSSD's
# secondary risk = set of cum LSSD's
# risk = the mLSSD
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
print(f"risk {risk} evaluation test {np.dot(prim, coef)}")

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
opt = {'title': "LSSD Port", 'tangent': True}
print("\n rate of returns vs risk representation")
rft = cr1.viewFrontiers(musharpe=0., randomport=100, options=opt)
print("\n sharpe vs rate of returns representation")
rft2 = cr1.viewFrontiers(data=rft, fig_type='Sharpe_RR')

#=============================================================================
# Sharpe vs. Sharpe2
# first Sharpe (default rtype)
cr1 = az.LSSDAnalyzer(coef, mktdata)
ww1 = cr1.getWeights(mu=0.)
RR1 = cr1.RR
risk1 = cr1.risk
prim1 = cr1.primary_risk_comp.copy()
seco1 = cr1.secondary_risk_comp.copy()
sharpe1 = cr1.sharpe
# second Sharpe2
cr2 = az.LSSDAnalyzer(coef, mktdata)
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
cr1 = az.LSSDAnalyzer(coef, mktdata)
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
cr1 = az.LSSDAnalyzer(coef, mktdata)
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
cr1 = az.LSSDAnalyzer(coef, mktdata)
ww1 = cr1.getWeights(mu=0.)
sharpe = cr1.sharpe
risk = cr1.risk

# compute RiskAverse portfolio for Lambda=sharpe
Lambda = sharpe
cr2 = az.LSSDAnalyzer(coef, mktdata)
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
#     crrx = az.LSSDAnalyzer(coef, mktdata, method=method)
#     toc = time.perf_counter()
#     wwx = crrx.getWeights(mu=0.)
#     tic = time.perf_counter() - toc
#     print(f"method: {method} time: {tic}")
#     xta[method] = pd.Series([tic], index=["Time"]).append(wwx)

# res = pd.DataFrame(xta)
# print(res.round(4))

#=============================================================================
# Example of rebalancing positions
cr1 = az.LSSDAnalyzer(coef, mktdata)

# existing positions and cash
ns = pd.Series(100, index=symb)
cash = 0.

# new positions and rolling info
pos = cr1.getPositions(mu=0., rtype='Sharpe', nshares=ns, cash=0.)
print(f" New position report\n {pos}")
