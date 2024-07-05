import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# create the SDAnalyzer class calculator
calc = az.MVAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='MinRisk' -> minimum SD portfolio
ww = calc.getWeights(rtype='MinRisk')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio SD
SD = calc.risk