# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 15:23:34 2021

@author: mirce
"""
# Examples
import numpy as np
import pandas as pd

import azapy as az

#=============================================================================
# Collect some market data
sdate = pd.to_datetime("2012-01-01")
edate = pd.to_datetime('today')
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdir = "./MkTdata"

# force=True read from alphavantage server
# force=False read from local directory if data exists
mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                     dir=mktdir, force=False) 

#=============================================================================
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
# risk = weighted sum of LSSD's 
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
# the same rate of returns.
test_ww1 = cr1.getWeights(mu=RR, rtype='Risk')
ww_comp = pd.DataFrame({"ww1": ww1, "test_ww1": test_ww1,
                        'diff': ww1-test_ww1})
print(f"Test for weights computation\n {ww_comp}")

#=============================================================================
#Frontier evaluations
print("\nFrontiers evaluations\n")
opt = {'title': "Portfolio frontiers", 'tangent': True}
print("\n rate of returns vs risk representation")
saveto = "graphics/frontiers_1.png"
rft = cr1.viewFrontiers(musharpe=0., randomport=100, options=opt, saveto=saveto)
print("\n sharpe vs rate of returns representation")
saveto = "graphics/frontiers_2.png"
rft2 = cr1.viewFrontiers(data=rft, fig_type='Sharpe_RR', saveto=saveto)