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
# rtype='Sharpe2' -> minimization of inverse mMAD-Sharpe ratio
# risk free rate set to 0 (default mu0=0)
ww = calc.getWeights(rtype='Sharpe2')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mMAD-Sharpe ratio
sharpe = calc.sharpe
# portfolio mMAD
mMAD = calc.risk
# delta-risk components
delta = calc.primary_risk_comp.copy()