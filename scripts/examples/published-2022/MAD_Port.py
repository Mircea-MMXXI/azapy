import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
# set the start and end dates of the time series
sdate = '2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, sdate=sdate, edate=edate)

# define the parameters of mMAD measure
# coefficients - the normalization is done internally
# implicitly the mMAD level is len(coef)
coef = [1] * 3

# create Port_MAD calculator (quarterly rebalanced)
# pname="mMAD-Sharpe Port" - portfolio name
calc = az.Port_MAD(mktdata, pname="mMAD-Sharpe Port")

# set the model parameters:
#    coef - m-level coefficients
#    mMAD-Sharpe optimal portfolio (default rtype='Sharpe')
#    risk free rate set to 0 (default mu0=0)
# returns the portfolio time series (backtesting)
port = calc.set_model(coef=coef)

# plot the portfolio time series + EMA30 and EMA200
fig = '../figs/Port_mMAD_Sharpe.png'
_ = calc.port_view(saveto=fig)
