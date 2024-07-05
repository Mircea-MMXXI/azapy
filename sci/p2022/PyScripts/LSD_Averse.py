import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mLSD measure
# coefficients - the normalization is done internally
# implicitly the mLSD level is len(coef)
coef = [1] * 3

# create the LSDAnalyzer class calculator
calc = az.LSDAnalyzer(coef, mktdata)

# compute optimal portfolio weights
# aversion=2 -> risk-aversion factor set to 2
# rtype='RiskAverse' -> optimal for fixed risk-aversion
ww = calc.getWeights(rtype='RiskAverse', aversion=2)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mLSD
mLSD = calc.risk
# delta-risk components
delta = calc.primary_risk_comp.copy()