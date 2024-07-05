import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mMAD measure
# coefficients - the normalization is done internally
# implicitly the mMAD level is len(coef)
coef = [1] * 3

# create the MADAnalyzer class calculator
calc = az.MADAnalyzer(coef, mktdata)

# compute optimal portfolio weights
# rtype='Risk' -> minimization of risk
# mu=0.035 -> portfolio expected rate of return set to 3.5%
ww = calc.getWeights(rtype='Risk', mu=0.035)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mMAD
mMAD = calc.risk
# delta-risk components
delta = calc.primary_risk_comp.copy()