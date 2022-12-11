import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
# set the end date of the time series
# default start date is sdate='2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, edate=edate)

# define the parameters of the mixture CVaR measure
# confidence levels
alpha = [0.975, 0.95, 0.9]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)

# create the CVaRAnalyzer class calculator
calc = az.CVaRAnalyzer(alpha, coef, mktdata)

# plot expected rate of return vs mCVaR
# set some plot features
opt = {'title': "mCVaR Port",
       'xlabel': "mCVaR : " + str(alpha) + ", equal weigthed",
       'tangent': True}
fig1 = '../figs/mCVaRFrontier1.png'
rft = calc.viewFrontiers(randomport=100, options=opt,
                         saveto=fig1)

# plot mCVaR-Sharpe vs expected rate of return
# set some plot features
opt2 = {'title': "mCVaR Port",
        'ylabel': "mCVaR-Sharpe"}
fig2 = '../figs/mCVaRFrontier2.png'
# reuse the previous plot computations (data=rft)
_ = calc.viewFrontiers(data=rft, fig_type='Sharpe_RR',
                       options=opt2, saveto=fig2)
