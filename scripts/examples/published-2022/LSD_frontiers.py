import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
# set the end date of the time series
# default start date is sdate='2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, edate=edate)

# define the parameters of mLSD measure
# coefficients - the normalization is done internally
# implicitly the mMAD level is len(coef)
coef = [1] * 3

# create the LSDAnalyzer class calculator
calc = az.LSDAnalyzer(coef, mktdata)

# plot expected rate of return vs mLSD
# set some plot features
opt = {'title': "mLSD Port",
       'xlabel': "mLSD : [L=3, equal weighted]",
       'tangent': True}
fig1 = '../figs/mLSDFrontier1.png'
rft = calc.viewFrontiers(musharpe=0, randomport=100,
                         options=opt, saveto=fig1)
                         
# plot mLSD-Sharpe vs expected rate of return
# set some plot features
opt2 = {'title': "mLSD Port",
        'ylabel': "mLSD-Sharpe"}
fig2 = '../figs/mLSDFrontier2.png'
# reuse the previous plot computations (data=rft)
_ = calc.viewFrontiers(data=rft, fig_type='Sharpe_RR',
                       options=opt2, saveto=fig2)
