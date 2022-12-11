import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mMAD measure
# coefficients - the normalization is done internally
# implicitly the mMAD level is len(coef)
coef = [1] * 3

# create the MADAnalyzer class calculator
calc = az.MADAnalyzer(coef, mktdata)

# compute optimal portfolio weights
# aversion=5 -> risk-aversion set to 5
# rtype='RiskAverse' -> optimal for fixed risk-aversion
ww = calc.getWeights(rtype='RiskAverse', aversion=5)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mMAD
mMAD = calc.risk
# delta-risk components
delta = calc.primary_risk_comp.copy()