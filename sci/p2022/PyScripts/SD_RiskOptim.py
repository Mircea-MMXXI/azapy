import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# create the SDAnalyzer class calculator
calc = az.SDAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='Risk' -> minimization of risk
# mu=0.035 -> portfolio expected rate of return set to 3.5%
ww = calc.getWeights(rtype='Risk', mu=0.035)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio SD
SD = calc.risk