# Compare BTSD for alpha0=0 with LSSD first order
import numpy as np
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Set the BTSD parameter alpha0
alpha0 = 0.00

#=============================================================================
# Compute Sharpe optimal portfolio
# build the analyzer object
cr1 = az.BTSDAnalyzer(alpha0, mktdata)
# computes Sharpe weights for 0 risk-free rate
ww1 = cr1.getWeights(mu=0.)
# print portfolio characteristics
# primary risk = [Delta-risk] (redundant)
# secondary risk = [Delta-risk] (redundant)
# risk = Delta-risk
# Share = Omega ratio
RR1 = cr1.RR
risk1 = cr1.risk
prim1 = cr1.primary_risk_comp.copy()
seco1 = cr1.secondary_risk_comp.copy()
sharpe1 = cr1.sharpe
print("\nSharpe optimal portfolio\n")
print(f"status {cr1.status}")
print(f"coef {ww1}")
print(f"Secondary risk {seco1}")
print(f"Primary risk {prim1}")
print(f"Sharpe {sharpe1}")
print(f"RR {RR1}")
print(f"risk {risk1}")

#=============================================================================
# Define mLSSD measure parameters coef
coef = [1]

#=============================================================================
# Compute Sharpe optimal portfolio
# build the analyzer object
cr2 = az.LSSDAnalyzer(coef, mktdata)
# computes Sharpe weights for 0 risk-free rate
ww2 = cr2.getWeights(mu=0.)
# print portfolio characteristics
# primary risk = set of MAD's
# secondary risk = set of cumulative MAD's
# risk = mMAD value
RR2 = cr2.RR
risk2 = cr2.risk
prim2 = cr2.primary_risk_comp.copy()
seco2 = cr2.secondary_risk_comp.copy()
sharpe2 = cr2.sharpe
print("\nSharpe optimal portfolio\n")
print(f"status {cr2.status}")
print(f"coef {ww2}")
print(f"Secondary risk {seco2}")
print(f"Primary risk {prim2}")
print(f"Sharpe {sharpe2}")
print(f"RR {RR2}")
print(f"risk {risk2} evaluation test {np.dot(prim2, coef)}")

#=============================================================================
# comparison
print(f"risk: BTSD {risk1} LSSD {risk2} diff {risk1-risk2}")