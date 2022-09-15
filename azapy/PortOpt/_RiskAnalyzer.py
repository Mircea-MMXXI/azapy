import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import warnings

from azapy.MkT.MkTcalendar import NYSEgen

_WW_TOL = 1.e-10

class _RiskAnalyzer:
    """
    Base class for all XXXAnalyzer classes. \n
    All derive classes need to implement: \n
        _risk_calc(self, prate, alpha) \n
        _risk_min(self, d=1) \n
        _sharpe_max(self) \n
        _sharpe_inv_min(self) \n
        _rr_max(self) \n
        _risk_averse(self) \n
        _risk_diversification \n

    Methods:
        * getWeights
        * getRisk
        * getPositions
        * getRiskComp
        * getDiversification
        * viewForntiers
        * set_rrate
        * set_mktdata
        * set_rtype
        * set_random_seed
        
    Public members:
        
    """

    def __init__(self, mktdata=None, colname='adjusted',
                 freq='Q', hlength=3.25, calendar=None, rtype='Sharpe'):
        """
        Constructor

        Parameters
        ----------
        `mktdata` : `pandas.DataFrame`, optional
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        `colname` : string, optional
            Name of the price column from mktdata used in the weights
            calibration. The default is 'adjusted'.
        `freq` : string, optional
            Rate of return horizon in number of business day. it could be
            'Q' for quarter or 'M' for month. The default is `'Q'`.
        `hlength` : float, optional
            History length in number of years used for calibration. A
            fractional number will be rounded to an integer number of months.
            The default is `3.25` year.
        `calendar` : `numpy.busdaycalendar`, optional
            Business days calendar. If is it `None` then the calendar will be 
            set to NYSE business calendar. The default is `None`.
        `rtype` : string, optional
            Optimization type. Possible values \n
                'Risk' : optimal portfolio for targeted expected rate of return.\n
                'Sharpe' : optimal Sharpe portfolio - maximization solution.\n
                'Sharpe2' : optimal Sharpe portfolio - minimization solution.\n
                'MinRisk' : minimum risk portfolio.\n
                'InvNrisk' : optimal portfolio with the same risk value as a 
                benchmark portfolio (e.g. equal weighted portfolio).\n
                'RiskAverse' : optimal portfolio for a fixed risk-aversion
                factor.\n
                `MaxDivers` : portfolio with maximum diversification.\n
                `Divers` : optimal diversified portfolio for targeted
                expected rate of return.\n
            The default is 'Sharpe'.

        Returns
        -------
        The object.
        """
        self.ww = None
        self.status = None

        self.last_data = None
        self.rrate = None
        self.nn = None
        self.mm = None
        self.muk = None
        self.mu = None
        
        self.alpha = [1.]
        self.coef = [1.]
        
        self.Lambda = None

        self.risk = None
        self.primary_risk_comp = None
        self.secondary_risk_comp = None
        self.sharpe = None
        self.RR = None

        self.rtype = None

        self.rng = None
       
        self.method = None
        
        self.divers = None
        self.risk_comp = None
        self._flag_risk_comp_calc = False
        
        self.set_rtype(rtype)

        if mktdata is not None:
            self.set_mktdata(mktdata, colname, freq, hlength)
        self.set_random_seed()

        
    def getWeights(self, rtype=None, mu=None, d=1, mu0=0., aversion=None,
                   ww0=None, rrate=None ):  
        """
        Computes the optimal portfolio weights.

        Parameters
        ----------
        `rtype` : `str`, optional;
            Optimization type. If is not `None` it will overwrite the value
            set by the constructor. The default is `None`.
        `mu` : float, optional;
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype` is equalt to `'Risk'` or `'Divers'`.
            The default is `None`.
        `d` : int, optional;
            Frontier type. Active only if `rtype='Risk'`. A value of `1` will
            trigger the evaluation of optimal portfolio along the efficient
            frontier. Otherwise, it will find the portfolio with the lowest
            rate of return along the inefficient portfolio frontier.
            The default is `1`.
        `mu0` : float, optional;
            Risk-free rate accessible to the investor.
            Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        `aversion` : float, optional;
            The value of the risk-aversion coefficient.
            Must be positive. Relevant only if `rtype='RiskAverse'`.
            The default is `None`.
        `ww0` : `list`, `numpy.array` or `pandas.Series`, optional;
            Targeted portfolio weights. 
            Relevant only if `rype='InvNrisk'`.
            Its length must be equal to the number of symbols in `rrate` (mktdata). 
            All weights must be >= 0 with sum > 0.
            If it is a `list` or a `numpy.array` then the weights are assumed to
            by in order of `rrate.columns`. If it is a `pandas.Series` then 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (same symbols, not necessary in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.
        `rrate` : `pandas.DataFrame`, optional;
            The portfolio components historical rates of returns.
            If it is not `None`, it will overwrite the `rrate` computed in the
            constructor from mktdata. The default is `None`. 

       Returns
        -------
        `pandas.Series` with portfolio weights.
        """
        self.ww = None
        self.risk = None
        self.primary_risk_comp = None
        self.secondary_risk_comp = None
        self.sharpe = None
        self.RR = None
        
        if rrate is not None:
            self.set_rrate(rrate)

        if not rtype is None:
            self.set_rtype(rtype)
            
        if self.rtype == "Risk":
            if mu is None:
                raise ValueError("for rtype='Risk' mu must have a value")
            elif mu > self.muk.max():
                self.mu = self.muk.max()
            elif mu < self.muk.min():
                self.mu = self.muk.min()
            else:
                self.mu = mu

            self._risk_min(d=d)
            
        elif self.rtype == 'Sharpe':           
            if mu0 is None:
                raise ValueError("for rtype='Sharpe' mu0 must have a value")
            if mu0 > self.muk.max():
                self.mu = self.muk.max() * 0.999
            else:
                self.mu = mu0

            self._sharpe_max()
            
        elif self.rtype == 'Sharpe2':          
            if mu0 is None:
                raise ValueError("for rtype='Sharpe2' mu0 must have a value")
            if mu0 > self.muk.max():
                self.mu = self.muk.max() * 0.999
            else:
                self.mu = mu0

            self._sharpe_inv_min() 
            
        elif self.rtype == "MinRisk":
            self.mu = self.muk.min()
            self._risk_min(d=1)
            
        elif self.rtype == "InvNrisk":    
            if ww0 is None:
                ww = np.full(len(self.rrate.columns), 
                             1/len(self.rrate.columns))
            else:
                if isinstance(ww0, pd.core.series.Series):
                    ww = np.array(ww0[self.rrate.columns])
                else:
                    ww = np.array(ww0)
                if any(ww < 0.):
                    raise ValueError("All ww0 must be non-negative")
                sww = ww.sum()
                if sww <= 0:
                    raise ValueError("At least one ww0 must be non zero")
                ww = ww / sww
                
            self.getRisk(ww)
            self._rr_max()  
            
        elif self.rtype == "RiskAverse":          
            if aversion is None:
                raise ValueError("for rtype='RiskAverse'"
                                 + " aversion must have a value")
            elif aversion < 0:
                raise ValueError("aversion must be positive")
                
            self.Lambda = aversion
            self._risk_averse()
            
        elif self.rtype == "Divers":
            if mu is None:
                raise ValueError("for rtype='Divers' mu must have a value")
            elif mu > self.muk.max():
                self.mu = self.muk.max()
            elif mu < self.muk.min():
                self.mu = self.muk.min()
            else:
                self.mu = mu

            self.getRiskComp()
            self._risk_diversification(d=d)
            
        elif self.rtype == "MaxDivers":  
            self.mu = self.muk.min()
            self.getRiskComp()
            self._risk_diversification()
            
            
        else:
            raise ValueError("you should not be here")
        
        if self.status != 0:
            warnings.warn(f"Warning: status {self.status} for {self.rtype}")
            
        self.ww = pd.Series(self.ww, index=self.rrate.columns, dtype='float64')
        self._norm_weights()
        return self.ww


    def getRisk(self, ww, rrate=None):
        """
        Returns the value of the dispersion (risk) measure for a give 
        portfolio.

        Parameters
        ----------
        `ww` : `list`, `numpy.array` or `pandas.Series`;
            Portfolio weights. Its length must be equal to the number of
            symbols in `rrate` (mktdata). All weights must be >=0 with 
            sum > 0.
            If it is a `list` or a `numpy.array` then the weights are assumed to
            by in order of `rrate.columns`. If it is a `pandas.Series` than 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (not necessary in the same order).
        `rrate` : `pandas.DataFrame`, optional;
            Contains the portfolio components historical
            rates of returns. If it is not `None`, it will overwrite the
            `rrate` computed in the constructor from `mktdata`.
            The default is `None`.

        Returns
        -------
        float
            The dispersion (risk) measure value.
        """
        if rrate is not None:
            self.set_rrate(rrate)

        if isinstance(ww, pd.core.series.Series):
            ww = ww[self.rrate.columns]
        w = np.array(ww)
        if any(w < 0.):
            raise ValueError("All ww must be non-negative")
        sw = w.sum()
        if sw <= 0:
            raise ValueError("At least one ww must be positive")
        w = w / sw

        prate = np.dot(self.rrate, w)

        self.primary_risk_comp = []
        self.secondary_risk_comp = []
        for alpha in self.alpha:
            self.status, secondary, primary = self._risk_calc(prate, alpha)
            if self.status != 0:
                warnings.warn(f"status: {self.status}, \
                              wrong risk calc for alpha {alpha}")
                return np.nan

            self.primary_risk_comp.append(primary)
            self.secondary_risk_comp.append(secondary)
        self.risk = np.dot(self.primary_risk_comp, self.coef)
        self.primary_risk_comp = np.array(self.primary_risk_comp)
        self.secondary_risk_comp = np.array(self.secondary_risk_comp)
        self.RR = np.dot(w, self.muk)
        self.ww = w

        return self.risk

        
    def getPositions(self, nshares=None, cash=0, ww=None, rtype=None, 
                     mu=None, mu0=0., aversion=None, ww0=None, ):
        """
        Computes the rebalanced number of shares.

        Parameters
        ----------
        `nshares` : `panda.Series`, optional;
            Initial number of shares per portfolio component.
            A missing component
            entry will be considered 0. A `None` value assumes that all
            components entries are 0. The name of the components must be
            present in the mrkdata. The default is `None`.
        `cash` : float, optional;
            Additional cash to be added to the capital. A
            negative entry assumes a reduction in the total capital
            available for rebalance. The total capital cannot be < 0.
            The default is 0. 
        `ww` : `panda.Series`, optional;
            External overwrite portfolio weights. 
            If it not set to `None` these
            weights will overwrite the calibrated.
            The default is `None`. 
        `rtype` : str, optional;
            Optimization type. If is not `None` it will overwrite the value
            set by the constructor. The default is `None`.
        `mu` : float, optional
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype='Risk'`
            The default is `None`.
        `mu0` : float, optional;
            Risk-free rate accessible to the investor.
            Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        `aversion` : float, optional;
            The value of the risk-aversion coefficient.
            Must be positive. Relevant only if `rtype='RiskAvers'`.
            The default is `None`.
        `ww0` : list (also `np.array` or `pandas.Series`), optional;
            Targeted portfolio weights 
            Relevant only if `rype='InvNrisk'`.
            Its length must be equal to the number of
            symbols in `rrate` (mktdata). 
            All weights must be >= 0 with sum > 0.
            If it is a list or a `numpy.array` then the weights are assumed to
            by in order of `rrate.columns`. If it is a `pandas.Series` then the 
            index should be compatible with the `rrate.columns` or mktdata 
            symbols (same symbols, not necessary in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.

        Returns
        -------
        `pandas.DataFrame`: the rolling information.

        Columns:

            - "old_nsh" :
                initial number of shares per portfolio component and
                the additional cash. These are input values.
            - "new_nsh" :
                the new number of shares per component plus the residual
                cash (due to the rounding to an integer number of shares).
                A negative entry means that the investor needs to add more
                cash to cover for the roundup shortfall.
                It has a small value.
            - "diff_nsh" :
                the number of shares that needs to be both/sold in order
                to rebalance the portfolio positions.
            - "weights" :
                portfolio weights used for rebalancing. The cash entry is
                the new portfolio value (invested capital).
            - "prices" :
                the share prices used for rebalances evaluations.

        Note: Since the prices are closing prices, the rebalance can be
        computed after the market close and before the 
        trading execution (next day). 
        Additional cash slippage may occur due
        to share price differential between the previous day closing and
        execution time.
        """
       
        if rtype is not None:
            self.set_rtype(rtype)

        ns = pd.Series(0, index=self.rrate.columns)
        if nshares is not None:
            ns = ns.add(nshares, fill_value=0)

        if len(self.rrate.columns) != len(ns):
            raise ValueError("nshares must by a subset "
                             f"of {self.rrate.columns}")

        pp = self.last_data[self.rrate.columns.to_list()]

        cap = ns.dot(pp) + cash
        if ww is None:
            ww = self.getWeights(rtype=rtype, mu=mu, mu0=mu0, ww0=ww0, 
                                 aversion=aversion)
        else:
            wwp = pd.Series(0, index=self.rrate.columns)
            ww = wwp.add(ww, fill_value=0)
            if len(self.rrate.columns) == len(ns):
                raise ValueError(f"ww must be among {self.rrate.columns}")

        newns = (ww / pp * cap).round(0)
        newcash = cap - newns.dot(pp)

        ns['cash'] = cash
        newns['cash'] = newcash
        ww['cash'] = cap - newcash
        pp['cash'] = 1.

        res = pd.DataFrame({'old_nsh': ns,
                            'new_nsh': newns,
                            'diff_nsh': newns - ns,
                            'weights': ww.round(6),
                            'prices': pp})

        return res


    def set_rrate(self, rrate):
        """
        Sets portfolio components historical rates of returns.
        It will overwrite the value computed by the constructor from mktdata.

        Parameters
        ----------
        `rrate` : `pandas.DataFrame`;
            Portfolio components historical rates of returns. The
            columns are: "date", "symbol1", "symbol2", etc.
        Returns
        -------
        None
        """
        self._flag_risk_comp_calc = False
        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean()
        self.rrate = rrate - self.muk


    def set_mktdata(self, mktdata, colname='adjusted',
                    freq='Q', hlength=3.25, calendar=None):
        """
        Sets historical market data. It will overwrite the choice made in the
        constructor.

        Parameters
        ----------
        `mktdata` : `pandas.DataFrame`;
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function.
        `colname` : `str`, optional
            Name of the price column from mktdata used in the weights
            calibration. The default is 'adjusted'.
        `freq` : `str`, optional;
            Rate of returns horizon in number of business day. It could be
            'Q' for quarter or 'M' for month. The default is 'Q'.
        `hlength` : `float`, optional;
            History length in number of years used for calibration. A
            fractional number will be rounded to an integer number of months.
            The default is `3.25`.
        `calendar` : `numpy.busdaycalendar`, optional;
            Business days calendar. If is it `None` then the calendar will be 
            set to NYSE business calendar.
            The default is `None`.

        Returns
        -------
        None
        """
        if mktdata is None:
            return

        self.mktdata = mktdata

        periods = 63 if freq == 'Q' else 21

        if calendar is None:
            calendar = NYSEgen()

        edate = mktdata.index[-1]
        sdate = edate - pd.offsets.DateOffset(months=round(hlength * 12, 0))
        sdate = np.busday_offset(sdate.date(),
                                 0, roll='backward',  busdaycal=calendar)
        rrate = mktdata.pivot(columns='symbol', values=colname)[sdate:]
        self.last_data = rrate.iloc[-1]
        rrate = rrate.pct_change(periods=periods).dropna()

        self.set_rrate(rrate)


    def set_rtype(self, rtype):
        """
        Sets the optimization type. It will overwrite the value set in the
        constructor.

        Parameters
        ----------
        `rtype` : `str`
            Optimization type.
            
        Returns
        -------
        `None`
        """
        rtypes = ["Sharpe", "Risk", "MinRisk", "Sharpe2", "InvNrisk",
                  "RiskAverse", "Divers", "MaxDivers"]

        if not rtype in rtypes:
            raise ValueError(f"rtype must be one of {rtypes}")
        self.rtype = rtype


    def viewFrontiers(self, minrisk=True, efficient=20, inefficient=20, 
                      sharpe=True, musharpe=0.,
                      component=True, randomport=20, inverseN=True,
                      maxdivers=False, divers_efficient=0, divers_inefficient=0,
                      addport=None, fig_type='RR_risk', **opt):
        """
        Computes the elements of the portfolio frontiers.

        Parameters
        ----------
        `minrisk` : Boolean, optional;
            If it is `True` then the minimum risk portfolio will be visible. 
            The default is `True`.
        `efficient` : int, optional;
            Number of points along the optimal frontier (equally spaced along
            the rate of return axis). The default is `20`.
        `inefficient` : int, optional;
            Number of points along the inefficient frontier (equally spaced
            along the rate of returns axis). The default is `20`.
        `sharpe` : Boolean, optional;
            If it is `True` then the maximum Sharpe portfolio will be visible. 
            The default is `True`.
        `musharpe` : float, optional;
            Risk-free rate value used in the evaluation of
            generalized Sharpe ratio. The default is `0`.
        `component` : Boolean, optional;
            If it is `True` then the single component portfolios added to 
            the plot. The default is `True`.
        `randomport` : int, optional;
            Number of portfolios with random weights (inefficient) to be
            added to the plot (for reference). The default is `20`.
        `inverseN` : boolean, optional;
            If it is `True` then the equally weighted portfolio and the optimal 
            portfolio  with the same dispersion (risk) value are added to
            the plot. The default is `True`.
        `maxdivers`: Boolean, optional;
            If it is `True` then the maximum diversified portfolio will be 
            visible. The default is `True`.
        `divers_efficient`: int, optional;
            Number of points along the diversified efficient frontier
            (equally spaced along the rate of return axis).
            The default is 20.
        `divers_inefficient`: int, optional;
            Number of points along the diversified inefficient frontier
            (equally spaced along the rate of return axis).
            The default is 20.
        `addport` : `dict` or `pandas.DataFrame`, optional;
            The weights of additional portfolio to be added to the plot.
            If it is a `dict` then the keys are the labels, and the values are 
            list of weights in the order of `rrate.columns`. If it is a 
            `pandas.DataFrame` the index are the labels, and each row is a set 
            of weights. The columns names should match the symbols names. 
            The default is `None`.
        `fig_type` : `str`, optional;
            Graphical representation format. \n
                * `'RR_risk'` : expected rate of return vs risk
                * `'Sharpe_RR'` : sharpe vs expected rate of return
                * `'Divers_RR'` : diversification vs expected rate of return
                
            The default is `'RR_risk'`.
        `opt` : optional;
            Additonal parameters:\n
                * `'title'` : The default is 'Portfolio frontiers'
                * `'xlabel'` : The default is \n
                    - `'risk'` if `fig_type='RR_risk'`
                    - `'rate of returns'` otherwise
                * `'ylabel'` : The default is \n
                    - `'rate of returns'` if `fig_type='RR_risk'`
                    - `'sharpe'` if `fig_type='RR_sharpe'`
                    - `'diversification'` if `fig_type=RR_divers` 
                * `'tangent'` : Boolean flag. If set to `True` the tangent
                    (to sharpe point) is added. It has effect only  if
                    `fig_type='RR_risk'`. The default is `True`.
                * `saveto` : `str`, optional; \n
                    File name to save the figure. The extension dictates the format:
                    png, pdf, svg, etc. For more details see the `mathplotlib`
                    documentation for `savefig`. The default is `None`.
                * `data` : `dict`, optional; \n
                    Numerical data to construct the plot. If it is not `None` it
                    will take precedence and no other numerical evaluations will be
                    performed. It is meant to produce different plot representations
                    without reevaluations. The default is `None`.

        Returns
        -------
        `dict`
            Numerical data used to make the plots. It can be passed back to
            reconstruct the plots without reevaluations.
        """
        options = defaultdict(lambda: None)
        for k, v in opt.items():
            if v is None:
                next
            elif isinstance(v, dict) & (k != 'data'):
                options.update(v)
            else:
                options[k] = v
                
        saveto = options['saveto']
        del options['saveto']
        res = options['data']
        del options['data']
 
        # calc
        if res is None:
            res = defaultdict(lambda: None)
            
            minrisk_calc = False
            if minrisk | (efficient != 0) | (inefficient != 0):
                minrisk_calc = True
            
            maxdivers_calc = False
            self.getRiskComp()
            if maxdivers | (divers_efficient != 0) | (divers_inefficient !=0): 
                maxdivers_calc = True
  
            # min risk
            res['risk_min'] = defaultdict(lambda: None)
            res['risk_min']['risk_min'] = minrisk
            if minrisk_calc:
                ww_min = self.getWeights(rtype='MinRisk')
                risk_min = self.risk
                rr_min = self.RR
                res['risk_min']['risk'] = risk_min
                res['risk_min']['rr'] = rr_min
                res['risk_min']['ww'] = ww_min
                res['risk_min']['sharpe'] = (np.dot(ww_min, self.muk) - musharpe) / risk_min
                res['risk_min']['divers'] = 1 - risk_min / np.dot(ww_min, self.risk_comp)
                res['risk_min']['label'] = 'MinR'
    
            # efficient frontier
            res['efficient'] = defaultdict(lambda: None)
            res['efficient']['efficient'] = efficient
            if efficient > 0:
                rr_eff = np.linspace(rr_min, self.muk.max(), efficient)
                risk_eff = []
                ww_eff = []
                RR_eff = []
                sharpe_eff = []
                divers_eff = []
                for rr in rr_eff:
                    ww_eff.append(self.getWeights(rtype='Risk', mu=rr))
                    risk_eff.append(self.risk)
                    RR_eff.append(self.RR)
                    sharpe_eff.append((np.dot(self.ww, self.muk) - musharpe) / self.risk)
                    divers_eff.append(1 - self.risk / np.dot(self.ww, self.risk_comp))
                res['efficient']['risk'] = np.array(risk_eff)
                res['efficient']['rr'] = np.array(RR_eff)
                res['efficient']['ww'] = ww_eff
                res['efficient']['sharpe'] = np.array(sharpe_eff)
                res['efficient']['divers'] = np.array(divers_eff)
    
            # inefficient frontier
            res['inefficient'] = defaultdict(lambda: None)
            res['inefficient']['inefficient'] = inefficient
            if inefficient > 0 :
                rr_ineff = np.linspace(self.muk.min(), rr_min, inefficient)
                risk_ineff = []
                ww_ineff = []
                RR_ineff = []
                sharpe_ineff = []
                divers_ineff = []
                for rr in rr_ineff:
                    ww_ineff.append(self.getWeights(rtype='Risk', mu=rr, d=-1))
                    risk_ineff.append(self.risk)
                    RR_ineff.append(self.RR)
                    sharpe_ineff.append((np.dot(self.ww, self.muk) - musharpe) / self.risk)
                    divers_ineff.append(1 - self.risk / np.dot(self.ww, self.risk_comp))
                res['inefficient']['risk'] = np.array(risk_ineff)
                res['inefficient']['rr'] = np.array(RR_ineff)
                res['inefficient']['ww'] = ww_ineff
                res['inefficient']['sharpe'] = np.array(sharpe_ineff)
                res['inefficient']['divers'] = np.array(divers_ineff)
    
            # sharpe point
            res['sharpe'] = defaultdict(lambda: None)
            res['sharpe']['sharpe'] = sharpe
            if sharpe:
                ww_sharpe = self.getWeights(rtype='Sharpe', mu0=musharpe)
                rr_sharpe = self.RR
                risk_sharpe = self.risk
                sharpe_sharpe = self.sharpe
                divers_sharpe = 1 - risk_sharpe / np.dot(ww_sharpe, self.risk_comp)
                res['sharpe']['mu'] = musharpe
                res['sharpe']['risk'] = risk_sharpe
                res['sharpe']['rr'] = rr_sharpe
                res['sharpe']['ww'] = ww_sharpe
                res['sharpe']['sharpe'] = sharpe_sharpe
                res['sharpe']['divers'] = divers_sharpe
                res['sharpe']['label'] = 'Sharpe'
    
            # portfolio components
            res['component'] = defaultdict(lambda: None)
            res['component']['component'] = component
            if component:
                risk_comp = []
                rr_comp = []
                label_comp = []
                ww_comp = []
                sharpe_comp = []
                divers_comp = []
                for k in range(self.mm):
                    ww = [0.] * self.mm
                    ww[k] = 1.
                    risk_comp.append(self.getRisk(ww))
                    rr_comp.append(self.muk[k])
                    label_comp.append(self.muk.index[k])
                    ww_comp.append(ww)
                    sharpe_comp.append((self.muk[k] - musharpe) / self.risk)
                    divers_comp.append(1 - self.risk / self.risk_comp[k])
                res['component']['risk'] = np.array(risk_comp)
                res['component']['rr'] = np.array(rr_comp)
                res['component']['ww'] = ww_comp
                res['component']['label'] = label_comp
                res['component']['sharpe'] = np.array(sharpe_comp)
                res['component']['divers'] = np.array(divers_comp)
    
            # random portfolios
            res['randomport'] = defaultdict(lambda: None)
            res['randomport']['randomport'] = randomport
            if randomport > 0:
                risk_rp = []
                rr_rp = []
                ww_rp = []
                sharpe_rp = []
                divers_rp = []
                for _ in range(randomport):
                    ww = self._ww_gen()
                    risk_rp.append(self.getRisk(ww))
                    rr_rp.append(self.RR)
                    ww_rp.append(ww)
                    sharpe_rp.append((np.dot(self.ww, self.muk) - musharpe) / self.risk)
                    divers_rp.append(1 - self.risk / np.dot(self.ww, self.risk_comp))
                res['randomport']['risk'] = np.array(risk_rp)
                res['randomport']['rr'] = np.array(rr_rp)
                res['randomport']['ww'] = ww_rp
                res['randomport']['sharpe'] = np.array(sharpe_rp)
                res['randomport']['divers'] = np.array(divers_rp)
    
            # inverse N
            res['inverseN'] = defaultdict(lambda: None)
            res['inverseN']['inverseN'] = inverseN
            if inverseN:
                risk_n = []
                rr_n = []
                ww_n = []
                sharpe_n = []
                divers_n = []
                
                # 1/N portfolio
                wwn = np.full(self.mm, 1 / self.mm)
                risk_n.append(self.getRisk(wwn))
                rr_n.append(self.RR)
                ww_n.append(wwn)
                sharpe_n.append((np.dot(self.ww, self.muk) - musharpe) / self.risk)
                divers_n.append(1 - self.risk / np.dot(self.ww, self.risk_comp))
                
                # optimal 1/N risk portfolio
                ww_n.append(self.getWeights(rtype="InvNrisk"))
                risk_n.append(self.risk)
                rr_n.append(self.RR)
                sharpe_n.append((np.dot(self.ww, self.muk) - musharpe) / self.risk)
                divers_n.append(1 - self.risk / np.dot(self.ww, self.risk_comp))
                
                res['inverseN']['risk'] = np.array(risk_n)
                res['inverseN']['rr'] = np.array(rr_n)
                res['inverseN']['ww'] = ww_n
                res['inverseN']['label'] = ['1/N', 'InvNrisk']
                res['inverseN']['sharpe'] = np.array(sharpe_n)
                res['inverseN']['divers'] = np.array(divers_n)
                
            # MaxDivers
            res['maxdivers'] = defaultdict(lambda: None)
            res['maxdivers']['maxdivers'] = maxdivers
            if maxdivers_calc:
                ww_d = self.getWeights('MaxDivers')  
                res['maxdivers']['risk'] = self.risk 
                rr_maxdivers = self.RR 
                res['maxdivers']['rr'] = rr_maxdivers
                res['maxdivers']['ww'] = ww_d
                res['maxdivers']['divers'] = self.divers
                res['maxdivers']['label'] = 'MaxD'
                res['maxdivers']['sharpe'] = (np.dot(ww_d, self.muk) - musharpe) / self.risk
                
            # divers efficient frontier
            res['divers_efficient'] = defaultdict(lambda: None)
            res['divers_efficient']['divers_efficient'] = divers_efficient
            if divers_efficient > 0:
                rr_d_eff = np.linspace(rr_maxdivers, self.muk.max(), divers_efficient)
                risk_d_eff = []
                ww_d_eff = []
                RR_d_eff = []
                sharpe_d_eff = []
                divers_d_eff = []
                for rr in rr_d_eff:
                    ww_d_eff.append(self.getWeights(rtype='Divers', mu=rr))
                    risk_d_eff.append(self.risk)
                    RR_d_eff.append(self.RR)
                    sharpe_d_eff.append((np.dot(self.ww, self.muk) - musharpe) / self.risk)
                    divers_d_eff.append(1 - self.risk / np.dot(self.ww, self.risk_comp))
                res['divers_efficient']['risk'] = np.array(risk_d_eff)
                res['divers_efficient']['rr'] = np.array(RR_d_eff)
                res['divers_efficient']['ww'] = ww_d_eff
                res['divers_efficient']['sharpe'] = np.array(sharpe_d_eff)
                res['divers_efficient']['divers'] = np.array(divers_d_eff)
                
            # divers inefficient frontier
            res['divers_inefficient'] = defaultdict(lambda: None)
            res['divers_inefficient']['divers_inefficient'] = divers_inefficient
            if divers_inefficient > 0 :
                rr_d_ineff = np.linspace(self.muk.min(), rr_maxdivers, divers_inefficient)
                risk_d_ineff = []
                ww_d_ineff = []
                RR_d_ineff = []
                sharpe_d_ineff = []
                divers_d_ineff = []
                for rr in rr_d_ineff:
                    ww_d_ineff.append(self.getWeights(rtype='Divers', mu=rr, d=-1))
                    risk_d_ineff.append(self.risk)
                    RR_d_ineff.append(self.RR)
                    sharpe_d_ineff.append((np.dot(self.ww, self.muk) - musharpe) / self.risk)
                    divers_d_ineff.append(1 - self.risk / np.dot(self.ww, self.risk_comp))
                res['divers_inefficient']['risk'] = np.array(risk_d_ineff)
                res['divers_inefficient']['rr'] = np.array(RR_d_ineff)
                res['divers_inefficient']['ww'] = ww_d_ineff
                res['divers_inefficient']['sharpe'] = np.array(sharpe_d_ineff)
                res['divers_inefficient']['divers'] = np.array(divers_d_ineff)
                
            # addport
            res['addport'] = defaultdict(lambda: None)
            res['addport']['addport'] = False
            if addport is not None:
                if isinstance(addport, dict):
                    aport = pd.DataFrame().from_dict(addport, 'index', columns=self.rrate.columns)
                elif isinstance(addport, pd.core.frame.DataFrame):
                    aport = addport[self.rrate.columns]
                else:
                    raise ValueError("addport must be a dict or a pandas.DataFrame, see doc")
     
                def ff_addport(x):
                    risk_addport = self.getRisk(x)
                    RR_addport = self.RR
                    sharpe_addport = (np.dot(x, self.muk) - musharpe) / risk_addport
                    divers_addport = 1 - risk_addport / np.dot(x, self.risk_comp) 
                    return pd.Series([risk_addport, RR_addport, sharpe_addport, divers_addport],
                                     index=['risk', 'rr', 'sharpe', 'divers'])
                
                res_addport = aport.apply(lambda x: ff_addport(x), axis=1)
                
                res['addport']['addport'] = True
                res['addport']['risk'] = res_addport['risk'].to_list()
                res['addport']['rr'] = res_addport['rr'].to_list()
                res['addport']['sharpe'] = res_addport['sharpe'].to_list()
                res['addport']['divers'] = res_addport['divers'].to_list()
                res['addport']['ww'] = aport.to_numpy().tolist()
                res['addport']['label'] = aport.index.to_list()

        # Plot
        res['options'] = options
        res['saveto'] = saveto
        if fig_type == 'RR_risk':
            self._plot_f1(res)
        elif fig_type == 'Sharpe_RR':
            self._plot_f2(res)
        elif fig_type == 'Divers_RR':
            self._plot_f3(res)
        else:
            raise ValueError(f"fig_type={fig_type} ; must be 'RR_risk', 'Sharpe_RR' or 'Divers_RR'")
            
        return res


    def _plot_f1(self, res):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        _opt = {'title': "Portfolio frontiers", 'xlabel': 'risk',
                'ylabel': 'rate of return', 'tangent': True}
        if res['options'] is None:
            opt = defaultdict(lambda x=0: None, _opt)
        else:
            opt = defaultdict(lambda x=0: None, res['options'])
            for key in _opt.keys():
                if opt[key] is None: opt[key] = _opt[key]

        ax.set(title=opt['title'], xlabel=opt['xlabel'], ylabel=opt['ylabel'])

        if res['risk_min']['risk_min']:
            ax.scatter(x=res['risk_min']['risk'], y=res['risk_min']['rr'],
                       marker='D', c='g', s=36)

        if res['efficient']['efficient'] > 0:
            ax.plot(res['efficient']['risk'], res['efficient']['rr'],
                    color='lightblue', linewidth=2)

        if res['inefficient']['inefficient'] > 0:
            ax.plot(res['inefficient']['risk'], res['inefficient']['rr'],
                    color='red', linewidth=2)

        if res['component']['component']:
            lf = 1.01
            for k in range(len(res['component']['risk'])):
                ax.scatter(x=res['component']['risk'][k],
                           y=res['component']['rr'][k],
                           marker='s', c='b', s=16)
                ax.annotate(res['component']['label'][k],
                            (res['component']['risk'][k] * lf,
                             res['component']['rr'][k] * lf))

        if res['randomport']['randomport'] > 0:
            ax.scatter(x=res['randomport']['risk'], y=res['randomport']['rr'],
                       marker='o', s=9, edgecolors='b', alpha=0.5)

        if res['inverseN']['inverseN']:
            ax.scatter(x=res['inverseN']['risk'], y=res['inverseN']['rr'],
                      marker='x', c='g', s=36)
            lf = 1.01
            for k in range(len(res['inverseN']['label'])):
                ax.annotate(res['inverseN']['label'][k],
                            (res['inverseN']['risk'][k] * lf,
                             res['inverseN']['rr'][k] * lf))
                
        if res['maxdivers']['maxdivers']:
            ax.scatter(x=res['maxdivers']['risk'], y=res['maxdivers']['rr'],
                       marker='x', c='b', s=36)
            lf = 1.01
            ax.annotate(res['maxdivers']['label'],
                        (res['maxdivers']['risk'] * lf,
                         res['maxdivers']['rr'] * lf))
    
        if res['divers_efficient']['divers_efficient'] > 0:
            ax.plot(res['divers_efficient']['risk'], res['divers_efficient']['rr'],
                    color='lightblue', linewidth=1, linestyle=':')

        if res['divers_inefficient']['divers_inefficient'] > 0:
            ax.plot(res['divers_inefficient']['risk'], res['divers_inefficient']['rr'],
                    color='red', linewidth=1, linestyle=':')
            
        if res['addport']['addport']:
            lf = 1.01
            for k in range(len(res['addport']['risk'])):
                ax.scatter(x=res['addport']['risk'][k],
                           y=res['addport']['rr'][k],
                           marker='o', s=9, edgecolors='r', alpha=0.5)
                ax.annotate(res['addport']['label'][k],
                            (res['addport']['risk'][k] * lf,
                             res['addport']['rr'][k] * lf))

        xl1, xl2 = plt.xlim()
        yl1, yl2 = plt.ylim()
        lff = 1.02
        xl2 *= lff
        yl2 *= lff
        ax.set(xlim=(xl1, xl2), ylim=(yl1, yl2))

        if res['sharpe']['sharpe']:
            ax.scatter(x=res['sharpe']['risk'], y=res['sharpe']['rr'],
                       marker='D', c='g', s=26)
            if opt['tangent']:
                sharpe = res['sharpe']['sharpe']
                musharpe = res['sharpe']['mu']
                if np.abs(sharpe) < 1.e-5:
                    ax.axhline(res['sharpe']['rr'])
                else:
                    x1 = xl1
                    y1 = sharpe * x1 + musharpe
                    if y1 < yl1:
                        y1 = yl1
                        x1 = (y1 - musharpe) / sharpe
                    x2 = xl2
                    y2 = sharpe * x2 + musharpe
                    if y2 > yl2:
                        y2 = yl2
                        x2 = (y2 - musharpe) / sharpe
                    ax.plot([x1, x2], [y1, y2], color='black', linewidth=1)
                    
        if res['saveto'] is not None:
            plt.savefig(res['saveto'])
        plt.show()


    def _plot_f2(self, res):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        _opt = {'title': "Portfolio frontiers", 'xlabel': 'rate of return',
                'ylabel': 'sharpe'}
        if res['options'] is None:
            opt = defaultdict(lambda x=0: None, _opt)
        else:
            opt = defaultdict(lambda x=0: None, res['options'])
            for key in _opt.keys():
                if opt[key] is None: opt[key] = _opt[key]

        ax.set(title=opt['title'], xlabel=opt['xlabel'], ylabel=opt['ylabel'])

        if res['risk_min']['risk_min']:
            ax.scatter(x=res['risk_min']['rr'], y=res['risk_min']['sharpe'], marker='D', c='g', s=36)

        if res['efficient']['efficient'] > 0:
            ax.plot(res['efficient']['rr'], res['efficient']['sharpe'],
                    color='lightblue', linewidth=2)

        if res['inefficient']['inefficient'] > 0:
            ax.plot(res['inefficient']['rr'], res['inefficient']['sharpe'],
                    color='red', linewidth=2)

        if res['component']['component']:
            lf = 1.01
            for k in range(len(res['component']['risk'])):
                ax.scatter(x=res['component']['rr'][k], y=res['component']['sharpe'][k],
                           marker='s', c='b', s=16)
                ax.annotate(res['component']['label'][k],
                            (res['component']['rr'][k] * lf, res['component']['sharpe'][k] * lf))
            
        if res['sharpe']['sharpe']:
            ax.scatter(x=res['sharpe']['rr'], y=res['sharpe']['sharpe'],
                       marker='D', c='g', s=26)

        if res['randomport']['randomport'] > 0:
            ax.scatter(x=res['randomport']['rr'], y=res['randomport']['sharpe'],
                       marker='o', s=9, edgecolors='b', alpha=0.5 )

        if res['inverseN']['inverseN']:
            ax.scatter(x=res['inverseN']['rr'], y=res['inverseN']['sharpe'],
                       marker='x', c='g', s=36)
            lf = 1.01
            for k in range(len(res['inverseN']['label'])):
                ax.annotate(res['inverseN']['label'][k],
                            (res['inverseN']['rr'][k] * lf, res['inverseN']['sharpe'][k] * lf))
                
        if res['maxdivers']['maxdivers']:
            ax.scatter(x=res['maxdivers']['rr'], y=res['maxdivers']['sharpe'],
                       marker='x', c='b', s=36)
            lf = 1.01
            ax.annotate(res['maxdivers']['label'],
                       (res['maxdivers']['rr'] * lf, res['maxdivers']['sharpe'] * lf))

        if res['divers_efficient']['divers_efficient'] > 0:
            ax.plot(res['divers_efficient']['rr'], res['divers_efficient']['sharpe'],
                    color='lightblue', linewidth=1, linestyle=':')

        if res['divers_inefficient']['divers_inefficient'] > 0:
            ax.plot(res['divers_inefficient']['rr'], res['divers_inefficient']['sharpe'],
                    color='red', linewidth=1, linestyle=':')
            
        if res['addport']['addport']:
            lf = 1.01
            for k in range(len(res['addport']['risk'])):
                ax.scatter(x=res['addport']['rr'][k],
                           y=res['addport']['sharpe'][k],
                           marker='o', s=9, edgecolors='r', alpha=0.5)
                ax.annotate(res['addport']['label'][k],
                            (res['addport']['rr'][k] * lf,
                             res['addport']['sharpe'][k] * lf))
            
        xl1, xl2 = plt.xlim()
        yl1, yl2 = plt.ylim()
        lff = 1.02
        xl2 *= lff
        yl2 *= lff
        ax.set(xlim=(xl1, xl2), ylim=(yl1, yl2))

        if res['saveto'] is not None:
            plt.savefig(res['saveto'])
        plt.show()


    def _plot_f3(self, res):
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
        _opt = {'title': "Portfolio frontiers", 'xlabel': 'rate of return',
                'ylabel': 'diversification'}
        if res['options'] is None:
            opt = defaultdict(lambda x=0: None, _opt)
        else:
            opt = defaultdict(lambda x=0: None, res['options'])
            for key in _opt.keys():
                if opt[key] is None: opt[key] = _opt[key]
    
        ax.set(title=opt['title'], xlabel=opt['xlabel'], ylabel=opt['ylabel'])

        if res['risk_min']['risk_min']:
            ax.scatter(x=res['risk_min']['rr'], y=res['risk_min']['divers'], marker='D', c='g', s=36)
            lf = 0.95
            ax.annotate(res['risk_min']['label'],
                       (res['risk_min']['rr'] * lf, res['risk_min']['divers'] * lf))
    
        if res['efficient']['efficient'] > 0:
            ax.plot(res['efficient']['rr'], res['efficient']['divers'],
                    color='lightblue', linewidth=1, linestyle=':')
    
        if res['inefficient']['inefficient'] > 0:
            ax.plot(res['inefficient']['rr'], res['inefficient']['divers'],
                    color='red', linewidth=1, linestyle=':')
    
        if res['component']['component']:
            lf = 1.01
            for k in range(len(res['component']['risk'])):
                ax.scatter(x=res['component']['rr'][k], y=res['component']['divers'][k],
                           marker='s', c='b', s=16)
                ax.annotate(res['component']['label'][k],
                            (res['component']['rr'][k] * lf, res['component']['divers'][k] * lf))
    
        if res['sharpe']['sharpe']:
            ax.scatter(x=res['sharpe']['rr'], y=res['sharpe']['divers'],
                       marker='D', c='g', s=26)
            lf = 1.01
            ax.annotate(res['sharpe']['label'],
                       (res['sharpe']['rr'] * lf, res['sharpe']['divers'] * lf))
    
        if res['randomport']['randomport'] > 0:
            ax.scatter(x=res['randomport']['rr'], y=res['randomport']['divers'],
                       marker='o', s=9, edgecolors='b', alpha=0.5 )
    
        if res['inverseN']['inverseN']:
            ax.scatter(x=res['inverseN']['rr'], y=res['inverseN']['divers'],
                       marker='x', c='g', s=36)
            lf = 1.01
            for k in range(len(res['inverseN']['label'])):
                ax.annotate(res['inverseN']['label'][k],
                            (res['inverseN']['rr'][k] * lf, res['inverseN']['divers'][k] * lf))
                
        if res['maxdivers']['maxdivers']:
            ax.scatter(x=res['maxdivers']['rr'], y=res['maxdivers']['divers'],
                       marker='x', c='b', s=36)
            lf = 1.01
            ax.annotate(res['maxdivers']['label'],
                       (res['maxdivers']['rr'] * lf, res['maxdivers']['divers'] * lf))
            
        if res['divers_efficient']['divers_efficient'] > 0:
            ax.plot(res['divers_efficient']['rr'], res['divers_efficient']['divers'],
                    color='lightblue', linewidth=2)
    
        if res['divers_inefficient']['divers_inefficient'] > 0:
            ax.plot(res['divers_inefficient']['rr'], res['divers_inefficient']['divers'],
                    color='red', linewidth=2)
            
        if res['addport']['addport']:
            lf = 1.01
            for k in range(len(res['addport']['risk'])):
                ax.scatter(x=res['addport']['rr'][k],
                           y=res['addport']['divers'][k],
                           marker='o', s=9, edgecolors='r', alpha=0.5)
                ax.annotate(res['addport']['label'][k],
                            (res['addport']['rr'][k] * lf,
                             res['addport']['divers'][k] * lf))
            
        xl1, xl2 = plt.xlim()
        yl1, yl2 = plt.ylim()
        lff = 1.02
        xl2 *= lff
        yl2 *= lff
        ax.set(xlim=(xl1, xl2), ylim=(yl1, yl2))
    
        if res['saveto'] is not None:
            plt.savefig(res['saveto'])
        plt.show()


    def set_random_seed(self, seed = 42):
        """
        Sets the seed for Dirichlet random generator used in viewFrontiers
        to select the inefficient portfolios.

        Parameters
        ----------
        `seed` : int, optional;
            The random generator seed, in case you want to set it to a weird
            value other than 42 :). The default is `42`.

        Returns
        -------
        None
        """
        self.rng = np.random.RandomState(seed)


    def _ww_gen(self):
        return self.rng.dirichlet([0.5] * self.mm)
    
    
    def _set_lp_method(self, method):
        lp_methods = ['ecos', 'highs-ds', 'highs-ipm', 'highs',
                       'interior-point', 'glpk', 'cvxopt']
        if not method in lp_methods:
            raise ValueError("unknown method {method} - "
                             + f"must be one of {lp_methods}")
        self.method = method
        
        
    def _set_qp_method(self, method):
        qp_methods = ['ecos', 'cvxopt']
        if not method in qp_methods:
            raise ValueError("unknown method {method} - "
                             + f"must be one of {qp_methods}")
        self.method = method
        
        
    def _set_socp_method(self, method):
        socp_methods = ['ecos', 'cvxopt']
        if not method in socp_methods:
            raise ValueError("unknown method {method} - "
                             + f"must be one of {socp_methods}")
        self.method = method
        

    
    def getRiskComp(self):
        """
        Returns the risk of each portfolio component.

        Returns
        -------
        `pandas.Series`
        """
        if self._flag_risk_comp_calc: 
            return self.risk_comp
            
        srisk = np.zeros(self.rrate.shape[1])
        for i in range(len(srisk)):
            wws = np.zeros(len(srisk))
            wws[i] = 1
            srisk[i] = self.getRisk(wws)
        
        self.risk_comp = pd.Series(srisk, self.rrate.columns)
        self._flag_risk_comp_calc = True

        return self.risk_comp
        
        
    
    def getDiversification(self, ww, rrate=None):
        """
        Returns the value of the diversification factor for a give 
        portfolio.

        Parameters
        ----------
        `ww` : `list`, `numpy.array` or `pandas.Series`;
            Portfolio weights. Its length must be equal to the number of
            symbols in `rrate` (mktdata). All weights must be >=0 with 
            sum > 0.
            If it is a list or a `numpy.array` then the weights are assumed to
            be in order of `rrate.columns`. If it is a `pandas.Series` than 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (not necessary in the same order).
        `rrate` : `pandas.DataFrame`, optional;
            Contains the portfolio components historical
            rates of returns. If it is not `None` then it will overwrite the
            rrate computed in the constructor from mktdata.
            The default is `None`.

        Returns
        -------
        `float` : The diversification value. \n
        
        
        """
        if rrate is not None:
            self.set_rrate(rrate)

        if isinstance(ww, pd.core.series.Series):
            ww = ww[self.rrate.columns]
        w = np.array(ww)
        if any(w < 0.):
            raise ValueError("All ww must be non-negative")
        sw = w.sum()
        if sw <= 0:
            raise ValueError("At least one ww must be positive")
        w = w / sw
        
        self.getRiskComp()
            
        wrisk = self.getRisk(w)
        
        crisk = np.dot(self.risk_comp, w)

        if crisk == 0:
            self.divers = 1
        else:
            self.divers = 1 - wrisk / crisk
        
        return self.divers
    
    
    def _norm_weights(self):
        self.ww = self.ww.apply(lambda x: x if x >= _WW_TOL else 0)
        self.ww /= self.ww.sum()


    # to be implemented in the deriv class
    def _risk_calc(self, prate, alpha):
        pass
    def _risk_min(self, d=1):
        pass
    def _sharpe_max(self):
        pass
    def _sharpe_inv_min(self):
        pass
    def _rr_max(self):
        pass
    def _risk_averse(self):
        pass
    def _risk_diversification(self, d=1):
        pass
    