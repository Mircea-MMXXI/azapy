import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mLSD measure
# coefficients - the normalization is done internally
# implicitly the mLSD level is len(coef)
coef = [1] * 3

# create the LSDAnalyzer class calculator
calc = az.LSDAnalyzer(coef, mktdata)

# compute optimal portfolio weights
# rtype='Sharpe2' -> minimization of inverse mLSD-Sharpe ratio
# risk free rate set to 0 (default mu0=0)
ww = calc.getWeights(mu=0, rtype='Sharpe2')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mLSD-Sharpe ratio
sharpe = calc.sharpe
# portfolio mLSD
mLSD = calc.risk
# delta-risk components
delta = calc.primary_risk_comp.copy()