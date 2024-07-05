import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mixture CVaR measure:
# confidence levels
alpha = [0.99, 0.975, 0.95]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)

# create the CVaRAnalyzer class calculator
calc = az.CVaRAnalyzer(alpha, coef, mktdata)

# compute optimal portfolio weights
# rtype='Sharpe2' -> minimization of inverse mCVaR-Sharpe ratio
# risk free rate set to 0 (default mu0=0)
ww = calc.getWeights(rtype='Sharpe2')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mCVaR-Sharpe ratio
sharpe = calc.sharpe
# portfolio mixture CVaR
mCVaR = calc.risk
# CVaR components
CVaR = calc.primary_risk_comp.copy()
# associated VaR values
VaR = calc.secondary_risk_comp.copy()