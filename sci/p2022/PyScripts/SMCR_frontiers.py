import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
# set the end date of the time series
# default start date is sdate='2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, edate=edate)

# define the parameters of mixture SMCR measure
# confidence levels
alpha = [0.9, 0.875, 0.85]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)

# create the SMCRAnalyzer class calculator
calc = az.SMCRAnalyzer(alpha, coef, mktdata)

# plot expected rate of return vs mSMCR
# set some plot features
opt = {'title': "mSMCR Port",
       'xlabel': "mSMCR : " + str(alpha) + ", equal weighted",
       'tangent': True}
fig1 = '../figs/mSMCRFrontier1.png'
rft = calc.viewFrontiers(musharpe=0, randomport=100,
                         options=opt, saveto=fig1)

# plot mSMCR-Sharpe vs expected rate of return
# set some plot features
opt2 = {'title': "mSMCR Port",
        'ylabel': "mSMCR-Sharpe"}
fig2 = '../figs/mSMCRFrontier2.png'
# reuse the previous plot computations (data=rft)
rft2 = calc.viewFrontiers(data=rft, fig_type='Sharpe_RR',
                          options=opt2, saveto=fig2)
