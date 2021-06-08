# -*- coding: utf-8 -*-
"""
Created on Sun May 30 23:14:06 2021

@author: mircea
"""

    
# Examples
import numpy as np
import pandas as pd

import azapy as az

#=============================================================================
# Collect some market data
sdate = pd.Timestamp("2012-01-01").normalize()
edate = pd.Timestamp.today().normalize()
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']
#symb = ["IHI", "SPY", "VGT", "PSJ", "PGF"]

mktdir = "./scripts/analyzers/MkTdata"

# force=True read from alphavantage server
# force=False read from local directory if data exists
rprice = az.readMkT(symb, dstart = sdate, dend = edate, 
                    dir=mktdir, force=False) 

# prepare mkt data: compute the 3-month rolling rate of return for adjusted 
# prices (column) and rearrange in format "date", "symbol1", "symbol2", etc.
# set a shorter hist for these examples (computational time is proportional
# to the square of number of historical observation)
hsdate =  pd.Timestamp("2020-01-01").normalize()

rrate = rprice.loc[rprice.index >= hsdate] \
    .pivot(columns='symbol', values='adjusted') \
    .pct_change(periods=62).dropna()

#=============================================================================
# Compute Sharpe optimal portfolio
# build the analyzer object
cr1 = az.GINIAnalyzer(rrate)
# computes Sharpe weights for 0 risk-free rate
ww1 = cr1.getWeights(mu=0.)
# print portfolio characteristics
# primary risk = set of GINI
# secondary risk = set of GINI
# risk = weighted sum of GINI
RR = cr1.RR
risk = cr1.risk
prim = cr1.primery_risk_comp.copy()
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
print(f"Test for the risk computation {test_risk} = {risk}")

# Test the Sharpe weights by estimating an optimal portfolio with 
# the same rate of returns.
test_ww1 = cr1.getWeights(mu=RR, rtype='Risk')
ww_comp = pd.DataFrame({"Original": ww1, "Test": test_ww1})
print(f"Test for weights computation {ww_comp}")

#=============================================================================
#Frontier evaluations - may take some time
# print("\nFrontiers evaluations\n")
# opt ={'title': "New Port", 'tangent': True}
# file = 'fig1w.svg'
# print("\n rate of return vs CVaR representation")
# rft = cr1.viewFrontiers(musharpe=0, randomport=1, options=opt)
# print("\n Sharpe vs rate of return representation")
# rft2 = cr1.viewFrontiers(data=rft, fig_type='Sharpe_RR')

#=============================================================================
# Test Sharpe vs. Sharpe2
# first Sharpe (default rtype)
cr1 = az.GINIAnalyzer(rrate)
ww1 = cr1.getWeights(mu=0.)
pd.Series(ww1, index=rrate.columns)
RR1 = cr1.RR
risk1 = cr1.risk
prim1 = cr1.primery_risk_comp.copy()
seco1 = cr1.secondary_risk_comp.copy()
sharpe1 = cr1.sharpe
# second Sharpe2
cr2 = az.GINIAnalyzer(rrate)
ww2 = cr2.getWeights(mu=0., rtype="Sharpe2")
pd.Series(ww2, index=rrate.columns)
RR2 = cr2.RR
risk2 = cr2.risk
prim2 = cr2.primery_risk_comp.copy()
seco2 = cr2.secondary_risk_comp.copy()
sharpe2 = cr2.sharpe
# print comparison
print(f"status {cr2.status} = {cr1.status}")
ww_comp = pd.DataFrame({"ww2": ww2, "ww1": ww1})
print(f"coef {ww_comp}")
seco_comp = pd.DataFrame({"seco2": seco2, "seco1": seco1})
print(f"Secondary risk {seco_comp}")
prim_comp = pd.DataFrame({"prim2": prim2, "prim1": prim1})
print(f"Primary risk {prim_comp}")
print(f"RR {RR2} = {RR1}")
print(f"risk {risk2} = {risk1}")
print(f"Sharpe {sharpe2} = {sharpe1}")

# Speed of Sharpe vs Sharpe2 - may take some time
# %timeit cr2.getWeights(mu=mu, rtype='Sharpe')
# %timeit cr2.getWeights(mu=mu, rtype='Sharpe2')

#=============================================================================
# Test for InvNrisk
cr1 = az.GINIAnalyzer(rrate)
# compute the risk of a equaly weighted portfolio
ww = np.ones(len(symb))
ww = ww / np.sum(ww)
risk = cr1.getRisk(ww)
# compute the weights of InvNrisk
ww1 = cr1.getWeights(mu=0., rtype="InvNrisk")
RR1 = cr1.RR
# compute the optimal portfolio for RR1 targeted rate of return 
ww2 = cr1.getWeights(mu=RR1, rtype="Risk")
# print comparison results
print(f"risk: 1/N port {risk} = InvNrisk {cr1.risk}")
ww_comp = pd.DataFrame({"InvNrisk": ww1, "Optimal": ww2})
print(f"weights: InvNrisk = Optimal {ww_comp}")

#=============================================================================
# Test for MinRisk
cr1 = az.GINIAnalyzer(rrate)
# compute the MinRisk portfolio
ww1 = cr1.getWeights(mu=0., rtype="MinRisk")
# test
ww2 = cr1.getWeights(mu=0., rtype="Risk")
# print comparison 
ww_comp = pd.DataFrame({"MinRisk": ww1, "Test": ww2})
print(f"weights: MinRisk = Optimal {ww_comp}")

#=============================================================================
# # speed comparison for different LP methods
# # may take some time to complete 
# # you have to uncomment the lines below
# crx1 = az.GINIAnalyzer(rrate, method='highs-ds')
# wwx1 = crx1.getWeights(mu=0.)
# print(wwx1)
# crx2 = az.GINIAnalyzer(rrate, method='highs-ipm')
# wwx2 = crx2.getWeights(mu=0.)
# print(wwx2)
# crx3 = az.GINIAnalyzer(rrate, method='highs')
# wwx3 = crx3.getWeights(mu=0.)
# print(wwx3)
# crx4 = az.GINIAnalyzer(rrate, method='interior-point')
# wwx4 = crx4.getWeights(mu=0.)
# print(wwx4)

# %timeit crx1.getWeights(mu=0.)
# %timeit crx2.getWeights(mu=0.)
# %timeit crx3.getWeights(mu=0.)
# %timeit crx4.getWeights(mu=0.)
