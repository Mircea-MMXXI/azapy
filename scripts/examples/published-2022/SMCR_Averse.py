import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mixture SMCR measure
# confidence levels
alpha = [0.9, 0.875, 0.85]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)

# create the SMCRAnalyzer class calculator
calc = az.SMCRAnalyzer(alpha, coef, mktdata)

# compute optimal portfolio weights
# rtype='RiskAverse' -> optimal for fixed risk-aversion
# aversion=0.05 -> risk-aversion factor set to 0.05
ww = calc.getWeights(rtype='RiskAverse', aversion=0.05)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mixture SMCR
mSMCR = calc.risk
# SMCR components
SMCR = calc.primary_risk_comp.copy()
# associated SMVaR values
SMVaR = calc.secondary_risk_comp.copy()