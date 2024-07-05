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
# rtype='MinRisk' -> minimum mLSD portfolio
ww = calc.getWeights(rtype='MinRisk')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mLSD
mLSD = calc.risk
# delta-risk components
delta = calc.primary_risk_comp.copy()