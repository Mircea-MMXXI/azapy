import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
mktdata = az.readMkT(symb)

# create the GINIAnalyzer class calculator
calc = az.GINIAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='Sharpe' -> maximization of Gini-Sharpe ratio
# risk free rate set to 0 (default mu0=0)
ww = calc.getWeights(rtype='Sharpe')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio Gini-Sharpe ratio
sharpe = calc.sharpe
# Gini dispersion
Gini = calc.risk