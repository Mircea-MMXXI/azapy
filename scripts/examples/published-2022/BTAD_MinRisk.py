import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mixture BTAD
# threshold levels
alpha = [-0.01, 0, 0.01]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)

# create the BTADAnalyzer class calculator
calc = az.BTADAnalyzer(alpha, coef, mktdata)

# compute optimal portfolio weights
# rtype='MinRisk' -> minimum mBTAD portfolio
ww = calc.getWeights(rtype='MinRisk')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mBTAD
mBTAD = calc.risk
# BTAD components
BTAD = calc.primary_risk_comp.copy()
