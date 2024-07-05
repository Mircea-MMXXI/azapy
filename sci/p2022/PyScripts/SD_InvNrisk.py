import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# create the SDAnalyzer class calculator
calc = az.SDAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='InvNrisk' -> optimal portfolio
#   with same SD as equal weighted portfolio
ww = calc.getWeights(rtype='InvNrisk')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio SD
SD = calc.risk