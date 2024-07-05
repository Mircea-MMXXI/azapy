import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# create the GINIAnalyzer class calculator
calc = az.GINIAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='MinRisk' -> minimum Gini portfolio
ww = calc.getWeights(rtype='MinRisk')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# Gini dispersion
Gini = calc.risk