import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
# set the start and end dates of the time series
sdate = '2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, sdate=sdate, edate=edate)

# define the parameters of mixture CVaR measure:
# confidence levels
alpha = [0.975, 0.95, 0.9]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)

# create Port_CVaR calculator (quarterly rebalanced)
# pname="mCVaR-Sharpe Port" - portfolio name
calc = az.Port_CVaR(mktdata, pname="mCVaR-Sharpe Port")

# set the model parameters:
#    alpha - confidence levels
#    coef - mixture coefficients
#    mCVaR-Sharpe optimal portfolio (default rtype='Sharpe')
#    risk free rate set to 0 (default mu0=0)
# returns the portfolio time series (backtesting)
port = calc.set_model(alpha=alpha, coef=coef)

# plot the portfolio time series + EMA30 and EMA200
fig = '../figs/Port_mCVaR_Sharpe.png'
_ = calc.port_view(saveto=fig)
