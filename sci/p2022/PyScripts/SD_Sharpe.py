import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# create the SDAnalyzer class calculator
calc = az.SDAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='Sharpe' -> maximization of Sharpe ratio
# risk free rate set to 0 (default mu0=0)
ww = calc.getWeights(rtype='Sharpe')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio Sharpe ratio
sharpe = calc.sharpe
# portfolio SD
SD = calc.risk