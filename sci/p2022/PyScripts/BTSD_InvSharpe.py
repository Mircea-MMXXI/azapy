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
# rtype='Sharpe2' -> minimization of inverse Sortino ratio
# risk free rate set to 0 (default mu0=0)
ww = calc.getWeights(rtype='Sharpe2')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio Sortino ratio
Sortino = calc.sharpe
# portfolio mBTSD
mBTSD = calc.risk
# BTSD components
BTSD = calc.primary_risk_comp.copy()