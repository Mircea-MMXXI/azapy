import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mixture BTSD
# threshold levels
alpha = [-0.01, 0, 0.01]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)

# create the BTSDAnalyzer class calculator
calc = az.BTSDAnalyzer(alpha, coef, mktdata)

# compute optimal portfolio weights
# aversion=1 -> risk-aversion factor set to 1
# rtype='RiskAverse' -> optimal with fixed risk-aversion
ww = calc.getWeights(rtype='RiskAverse', aversion=1)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mBTSD
mBTSD = calc.risk
# BTSD components
BTSD = calc.primary_risk_comp.copy()
