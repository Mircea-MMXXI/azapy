import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
# set the end date of the time series
# default start date is sdate='2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, edate=edate)

# define the parameters of the mixture BTAD
# threshold levels
alpha = [-0.01, 0, 0.01]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)


# create the BTADAnalyzer class calculator
calc = az.BTADAnalyzer(alpha, coef, mktdata)

# plot expected rate of return vs mBTAD
# set some plot features
opt = {'title': "mBTAD Port",
       'xlabel': "mBTAD : " + str(alpha) + ", equal weighted",
       'tangent': True}
fig1 = '../figs/mBTADFrontier1.png'
rft = calc.viewFrontiers(musharpe=0, randomport=100,
                         options=opt, saveto=fig1)

# plot Omega vs expected rate of return
# set some plot features
opt2 = {'title': "mBTAD Port",
        'ylabel': "Omega"}
fig2 = '../figs/mBTADFrontier2.png'
# reuse the previous plot computations (data=rft)
_ = calc.viewFrontiers(data=rft, fig_type='Sharpe_RR',
                       options=opt2, saveto=fig2)
