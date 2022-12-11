import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import copy as copy
import warnings
import time

from azapy.MkT.MkTcalendar import NYSEgen

_WW_TOL_ = 1.e-10
_MU_ABS_TOL_ = 1.e-6
_RISK_REL_TOL_ = 0.999999
_RR_ABS_TOL_ = 1.e-4


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
        _risk_diversification(self, d=1) \n
        _risk_inv_diversification(self, d=1) \n
        _rr_max_diversification(self) \n
        _set_method(self, method) \n

    Methods:
        * getWeights
        * getRisk
        * getPositions
        * getRiskComp
        * getDiversification
        * viewForntiers
        * set_rrate
        * set_method
        * set_mktdata
        * set_rtype
        * set_random_seed
        
    Public members:
        
    """

    def __init__(self, mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, calendar=None, rtype='Sharpe', name=None):
        """
        Constructor

        Parameters
        ----------
        `mktdata` : `pandas.DataFrame`, optional;
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        `colname` : `str`, optional;
            Name of the price column from mktdata used in the weights
            calibration. The default is `'adjusted'`.
        `freq` : `str`, optional;
            Rate of return horizon in number of business day. It could be
            `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
        `hlength` : `float`, optional;
            History length in number of years used for calibration. A
            fractional number will be rounded to an integer number of months.
            The default is `3.25` years.
        `calendar` : `numpy.busdaycalendar`, optional;
            Business days calendar. If is it `None` then the calendar will be 
            set to NYSE business calendar. The default is `None`.
        `rtype` : `str`, optional;
            Optimization type. Possible values: \n
                `'Risk'` : optimal risk portfolio for targeted expected rate of 
                return.\n
                `'Sharpe'` : optimal Sharpe portfolio - maximization solution.\n
                `'Sharpe2'` : optimal Sharpe portfolio - minimization solution.\n
                `'MinRisk'` : minimum risk portfolio.\n
                `'RiskAverse'` : optimal risk portfolio for a fixed 
                risk-aversion factor.\n
                `'InvNrisk'` : optimal risk portfolio with the same risk value 
                as a benchmark portfolio (e.g. same as equal weighted 
                portfolio).\n
                `'Diverse'` : optimal diversified portfolio for targeted
                expected rate of return (max of inverse 1-Diverse).\n
                `'Diverse2'` : optimal diversified portfolio for targeted
                expected rate of return (min of 1-Diverse).\n
                `'MaxDiverse'` : maximum diversified portfolio.\n
                `'InvNdiverse'` : optimal diversified portfolio with the same
                diversification factor as a benchmark portfolio 
                (e.g. same as equal weighted portfolio).\n
                `'InvNdrr'` : optimal diversified portfolio with the same 
                expected rate of return as a benchmark portfolio
                (e.g. same as equal weighted portfolio).\n
            The default is `'Sharpe'`.
        `name` : `str`, optional;
            Object name. Deafult value is `None`.

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
        self.diverse = None
        self.time_level1 = None
        self.time_level2 = None
        self.time_level3 = None

        self.rtype = None
        self.rng = None
        self.method = None
        
        self.methods = None
        self.rtypes = None
        self.name = None
        
        self.risk_comp = None
        self._flag_risk_comp_calc = False
        self._status_risk_comp_calc = None
        
        self.set_rtype(rtype)
        
        self.name = name

        if mktdata is not None:
            self.set_mktdata(mktdata, colname, freq, hlength)
        self.set_random_seed()
        
        self.verbose = False
        
        
    def _reset_output(self):
        self.status = None
        self.ww = None
        self.risk = None
        self.primary_risk_comp = None
        self.secondary_risk_comp = None
        self.sharpe = None
        self.RR = None
        self.diverse = None
        self.time_level1 = None
        self.time_level2 = None
        self.time_level3 = None

        
    def getWeights(self, rtype=None, mu=None, d=1, mu0=0., aversion=None,
                   ww0=None, rrate=None, verbose=False):  
        """
        Computes the optimal portfolio weights.

        Parameters
        ----------
        `rtype` : `str`, optional;
            Optimization type. If is not `None` it will overwrite the value
            set by the constructor. The default is `None`.
            Other possible values: \n
                `'Risk'` : optimal risk portfolio for targeted expected rate of 
                return.\n
                `'Sharpe'` : optimal Sharpe portfolio - maximization solution.\n
                `'Sharpe2'` : optimal Sharpe portfolio - minimization solution.\n
                `'MinRisk'` : minimum risk portfolio.\n
                `'RiskAverse'` : optimal risk portfolio for a fixed 
                risk-aversion factor.\n
                `'InvNrisk'` : optimal risk portfolio with the same risk value 
                as a benchmark portfolio (e.g. same as equal weighted 
                portfolio).\n
                `'Diverse'` : optimal diversified portfolio for targeted
                expected rate of return (max of inverse 1-Diverse).\n
                `'Diverse2'` : optimal diversified portfolio for targeted
                expected rate of return (min of 1-Diverse).\n
                `'MaxDiverse'` : maximum diversified portfolio.\n
                `'InvNdiverse'` : optimal diversified portfolio with the same
                diversification factor as a benchmark portfolio 
                (e.g. same as equal weighted portfolio).\n
                `'InvNdrr'` : optimal diversified portfolio with the same 
                expected rate of return as a benchmark portfolio
                (e.g. same as equal weighted portfolio).\n
        `mu` : `float`, optional;
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype` is equalt to `'Risk'` or `'Divers'`.
            The default is `None`.
        `d` : `int`, optional;
            Frontier type. Active only if `rtype='Risk'`. A value of `1` will
            trigger the evaluation of optimal portfolio along the efficient
            frontier. Otherwise, it will find the portfolio with the lowest
            rate of return along the inefficient portfolio frontier.
            The default is `1`.
        `mu0` : `float`, optional;
            Risk-free rate accessible to the investor.
            Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        `aversion` : `float`, optional;
            The value of the risk-aversion coefficient.
            Must be positive. Relevant only if `rtype='RiskAverse'`.
            The default is `None`.
        `ww0` : `list`, `numpy.array` or `pandas.Series`, optional;
            Targeted portfolio weights. 
            Relevant only if `rype='InvNrisk'`.
            Its length must be equal to the number of symbols in `rrate` 
            (mktdata). All weights must be >= 0 with sum > 0.
            If it is a `list` or a `numpy.array` then the weights are assumed 
            to be in order of `rrate.columns`. If it is a `pandas.Series` then 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (same symbols, not necessary in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.
        `rrate` : `pandas.DataFrame`, optional;
            The portfolio components historical rates of returns.
            If it is not `None`, it will overwrite the `rrate` computed in the
            constructor from mktdata. The default is `None`. 
        `verbose` : Boolean, optiona;
            If it set to `True` then it will print messages when the optimal
            portfolio degenerates to a single asset portfolio as a limited 
            case. 
            The default is `False`.

       Returns
        -------
        `pandas.Series` with portfolio weights.
        """
        toc = time.perf_counter()
        self._reset_output()
        self.verbose = verbose
 
        d = 1 if d == 1 else -1

        if rrate is not None:
            self.set_rrate(rrate)
            
        if not rtype is None:
            self.set_rtype(rtype)
            
        # calc - grand dispatch 
        if self.rtype == "Risk":
            self._calc_Risk(mu, d)
        elif self.rtype == 'Sharpe':           
            self._calc_Sharpe(mu0)
        elif self.rtype == 'Sharpe2':          
            self._calc_Sharpe2(mu0)
        elif self.rtype == "MinRisk":
            self._calc_MinRisk()
        elif self.rtype == "InvNrisk":    
            self._calc_InvNrisk(ww0)
        elif self.rtype == "RiskAverse":   
            self._calc_RiskAverse(aversion) 
        elif self.rtype == "Diverse2":
            self._calc_Diverse2(mu, d)
        elif self.rtype == "Diverse":
            self._calc_Diverse(mu, d)
        elif self.rtype == "MaxDiverse":  
            self._calc_MaxDiverse()
        elif self.rtype == "InvNdiverse":
            self._calc_InvNdiverse(ww0)
        elif self.rtype == "InvNdrr":
            self._calc_InvNdrr(ww0)
        else:
            raise ValueError("you should not be here - ever")

        if self.status != 0:
            warnings.warn(f"Warning: {self.rtype} on {self.rrate.index[-1]} "
                          f"status {self.status}")
        else:
            self.ww = pd.Series(self.ww, index=self.rrate.columns)
            self._norm_weights()
            
        self.time_level1 = time.perf_counter() - toc
        
        return self.ww
    
    
    def _calc_Risk(self, mu, d):
        mu_max = self.muk.max()
        mu_min = self.muk.min()
        if mu is None:
            raise ValueError("for rtype='Risk' mu must have a value")
        elif d == 1 and mu - mu_max > -np.abs(mu_max) * _MU_ABS_TOL_:
            self._single_asset_port('max')
            if self.verbose:
                print(f"rtype 'Risk' for mu {mu} for {self.name} on "
                      f"{self.rrate.index[-1]}, defaulted "
                      f"to single portfolio [{self.muk.idxmax()}]")
        elif d == -1 and mu - mu_min < np.abs(mu_min) * _MU_ABS_TOL_:
            self._single_asset_port('min')
            if self.verbose:
                print(f"rtype 'Risk' for mu {mu} for {self.name} on "
                      f"{self.rrate.index[-1]}, defaulted "
                      f"to single portfolio [{self.muk.idxmin()}]")
        else:
            self.mu = mu
            self._risk_min(d=d)
            
            
    def _calc_Sharpe(self, mu0):
        mu_max = self.muk.max()
        if mu0 is None:
            raise ValueError("for rtype='Sharpe' mu0 must have a value")
        if mu0 - mu_max >= -np.abs(mu_max) * _MU_ABS_TOL_:
            self._single_asset_port('max')
            if self.status == 0:
                self.sharpe = (self.RR - mu0) / self.risk
                if self.verbose:
                    print(f"rtype 'Sharpe' for mu0 {mu0} for {self.name} on "
                          f"{self.rrate.index[-1]}, defaulted "
                          f"to single portfolio [{self.muk.idxmax()}]")
        else:
            self.mu = mu0
            self._sharpe_max()
            
            
    def _calc_Sharpe2(self, mu0):
        mu_max = self.muk.max()
        if mu0 is None:
            raise ValueError("for rtype='Sharpe' mu0 must have a value")
        if mu0 - mu_max >= -np.abs(mu_max) * _MU_ABS_TOL_:
            self._single_asset_port('max')
            if self.status == 0:
                self.sharpe = (self.RR - mu0) / self.risk
                if self.verbose:
                    print(f"rtype 'Sharpe2' for mu0 {mu0} for {self.name} on "
                          f"{self.rrate.index[-1]}, defaulted "
                          f"to single portfolio [{self.muk.idxmax()}]")
        else:
            self.mu = mu0
            self._sharpe_inv_min() 
            
            
    def _calc_MinRisk(self):
        self._calc_Risk(self.muk.min(), d=1)
        
        
    def _calc_InvNrisk(self, ww0):
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
            
        self.getRiskComp()
        if self.status == 0:
            self.getRisk(ww)
            if self.status == 0:
                symb_max = self.muk.idxmax()
                if self.risk < self.risk_comp[symb_max] * _RISK_REL_TOL_:
                    self._rr_max()  
                else:
                    w = pd.Series(0., index=self.rrate.columns)
                    w[symb_max] = 1.
                    self.getRisk(w)
                    if self.verbose:
                        print(f"InvNrisk for {self.name} on "
                              f"{self.rrate.index[-1]}, defaulted "
                              f"to single portfolio [{symb_max}]")
                    
                    
    def _calc_RiskAverse(self, aversion):
        if aversion is None:
            raise ValueError("for rtype='RiskAverse'"
                             " aversion must have a positive value")
        elif aversion < 0:
            raise ValueError("aversion must be positive (>0)")
            
        self.Lambda = aversion
        self._risk_averse()
        
        
    def _calc_Diverse(self, mu, d):
        mu_max = self.muk.max()
        mu_min = self.muk.min()
        if mu is None:
            raise ValueError("for rtype='Diverse' mu must have a value")
        elif d == 1 and mu - mu_max > -np.abs(mu_max) * _MU_ABS_TOL_:
            self._single_asset_port('max')
            if self.verbose:
                print(f"rtype 'Diverse' for mu {mu} for {self.name} on "
                      f"{self.rrate.index[-1]}, defaulted "
                      f"to single portfolio [{self.muk.idxmax()}]")
        elif d == -1 and mu - mu_min < np.abs(mu_min) * _MU_ABS_TOL_:
            self._single_asset_port('min')
            if self.verbose:
                print(f"rtype 'Diverse' for mu {mu} for {self.name} on "
                      f"{self.rrate.index[-1]}, defaulted "
                      f"to single portfolio [{self.muk.idxmin()}]")
        else:
            self.mu = mu
            self.getRiskComp()
            if self.status == 0:
                self._risk_inv_diversification(d=d)
                
                
    def _calc_Diverse2(self, mu, d):
        mu_max = self.muk.max()
        mu_min = self.muk.min()
        if mu is None:
            raise ValueError("for rtype='Diverse2' mu must have a value")
        elif d == 1 and mu - mu_max > -np.abs(mu_max) * _MU_ABS_TOL_:
            self._single_asset_port('max')
            if self.verbose:
                print(f"rtype 'Diverse2' for mu {mu} for {self.name} on "
                      f"{self.rrate.index[-1]}, defaulted "
                      f"to single portfolio [{self.muk.idxmax()}]")
        elif d == -1 and mu - mu_min < np.abs(mu_min) * _MU_ABS_TOL_:
            self._single_asset_port('min')
            if self.verbose:
                print(f"rtype 'Diverse2' for mu {mu} for {self.name} on "
                      f"{self.rrate.index[-1]}, defaulted "
                      f"to single portfolio [{self.muk.idxmin()}]")
        else:
            self.mu = mu
            self.getRiskComp()
            if self.status == 0:
                self._risk_diversification(d=d)
                
                
    def _calc_MaxDiverse(self):
        self._calc_Diverse(self.muk.min(), d=1)
        
        
    def _calc_InvNdiverse(self, ww0):
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
            
        self.getDiversification(ww)
        if self.status == 0:
            self._rr_max_diversification()  
            
            
    def _calc_InvNdrr(self, ww0):
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
            
        mu = np.dot(ww, self.muk)
        self._calc_MaxDiverse()
        if self.status == 0:
            if mu > self.RR:
                self._calc_Diverse(mu, d=1)
            elif mu < self.RR:
                self._calc_Diverse(mu, d=-1)
        

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
            If it is a `list` or a `numpy.array` then the weights are assumed 
            to be in order of `rrate.columns`. If it is a `pandas.Series` than 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (not necessary in the same order).
        `rrate` : `pandas.DataFrame`, optional;
            Contains the portfolio components historical
            rates of returns. If it is not `None`, it will overwrite the
            `rrate` computed in the constructor from `mktdata`.
            The default is `None`.

        Returns
        -------
        `float`
            The dispersion (risk) measure value.
        """
        self._reset_output()
        toc = time.perf_counter()
        
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
        self.time_level1 = time.perf_counter() - toc
        
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
        `cash` : `float`, optional;
            Additional cash to be added to the capital. A
            negative entry assumes a reduction in the total capital
            available for rebalance. The total capital cannot be < 0.
            The default is 0. 
        `ww` : `panda.Series`, optional;
            External overwrite portfolio weights. 
            If it not set to `None` these
            weights will overwrite the calibration results.
            The default is `None`. 
        `rtype` : `str`, optional;
            Optimization type. If is not `None` it will overwrite the value
            set by the constructor. The default is `None`.
        `mu` : `float`, optional
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype='Risk'` and `rtype='Divers'`.
            The default is `None`.
        `mu0` : `float`, optional;
            Risk-free rate accessible to the investor.
            Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        `aversion` : `float`, optional;
            The value of the risk-aversion coefficient.
            Must be positive. Relevant only if `rtype='RiskAvers'`.
            The default is `None`.
        `ww0` : list (also `numpy.array` or `pandas.Series`), optional;
            Targeted portfolio weights 
            Relevant only if `rype='InvNrisk'`.
            Its length must be equal to the number of
            symbols in `rrate` (mktdata). 
            All weights must be >= 0 with sum > 0.
            If it is a `list` or a `numpy.array` then the weights are assumed to
            be in order of `rrate.columns`. If it is a `pandas.Series` then the 
            index should be compatible with the `rrate.columns` or mktdata 
            symbols (same symbols, not necessary in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.

        Returns
        -------
        `pandas.DataFrame`: the rolling information.

        Columns:

            - `"old_nsh"` :
                initial number of shares per portfolio component and
                the additional cash. These are input values.
            - `"new_nsh"` :
                the new number of shares per component plus the residual
                cash (due to the rounding to an integer number of shares).
                A negative entry means that the investor needs to add more
                cash to cover for the roundup shortfall.
                It has a small value.
            - `"diff_nsh"` :
                the number of shares that needs to be both/sold in order
                to rebalance the portfolio positions.
            - `"weights"` :
                portfolio weights used for rebalancing. The cash entry is
                the new portfolio value (invested capital).
            - `"prices"` :
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
            raise ValueError("nshares must be a subset " + \
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
        `colname` : `str`, optional;
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
        `None`
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
        `rtype` : `str`;
            Optimization type.
            
        Returns
        -------
        `None`
        """
        self.rtypes = ['Risk', 'MinRisk', 'Sharpe', 'Sharpe2', 'RiskAverse',
                       'InvNrisk', 'Diverse', 'Diverse2', 'MaxDiverse', 
                       'InvNdiverse', 'InvNdrr']

        if not rtype in self.rtypes:
            raise ValueError(f"rtype (value {rtype}) must be "
                             f"one of {self.rtypes}")
        self.rtype = rtype


    def viewFrontiers(self, minrisk=True, efficient=20, inefficient=20, 
                      sharpe=True, musharpe=0., maxdiverse=False, 
                      diverse_efficient=0, diverse_inefficient=0, invN=True, 
                      invNrisk=True, invNdiverse=False, invNdrr=False,
                      component=True, randomport=20, addport=None, 
                      fig_type='RR_risk', **opt):
        """
        Computes the elements of the portfolio frontiers.

        Parameters
        ----------
        `minrisk` : Boolean, optional;
            If it is `True` then the minimum risk portfolio will be visible. 
            The default is `True`.
        `efficient` : `int`, optional;
            Number of points along the efficient risk frontier 
            (equally spaced along the rate of return axis). 
            The default is `20`.
        `inefficient` : `int`, optional;
            Number of points along the inefficient risk frontier 
            (equally spaced along the rate of returns axis). 
            The default is `20`.
        `sharpe` : Boolean, optional;
            If it is `True` then the maximum Sharpe portfolio will be visible. 
            The default is `True`.
        `musharpe` : `float`, optional;
            Risk-free rate value used in the evaluation of
            generalized Sharpe ratio. The default is `0`.
        `maxdiverse`: Boolean, optional;
            If it is `True` then the maximum diversified portfolio will be 
            visible. The default is `True`.
        `diverse_efficient`: `int`, optional;
            Number of points along the efficient diversification frontier
            (equally spaced along the rate of return axis).
            The default is `20`.
        `diverse_inefficient`: `int`, optional;
            Number of points along the inefficient diversification frontier
            (equally spaced along the rate of return axis).
            The default is `20`.
        `invN` : Boolean, optional;
            If it is `True`, then the equal weighted portfolio and the optimal 
            portfolio with the same risk value are added to
            the plot. The default is `True`.
        `invNrisk` : Boolena, optional;
            If it is `True`, then the efficient risk portfolio with same risk 
            as equal weighth portfolio is added to the plot.
            The defualt is `False`.
        `invNdiverse` : Boolean, optional:
            If it is `True`, then the efficient diversified portfolio with the 
            same diversification factor as the equal weighted portfolio is
            added to the plot. The default is `False`.
        `invNdrr` : Boolena, optional;
            If it is `True`, then the efficient diversified portfolio with the 
            same expected rate of return as the equal weighted portfolio is
            added to the plot. The default value is `False`.
        `component` : Boolean, optional;
            If it is `True`, then the single component portfolios are added to 
            the plot. The default is `True`.
        `randomport` : `int`, optional;
            Number of portfolios with random weights (inefficient) to be
            added to the plot (for reference). The default is `20`.
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
                * `'Diverse_RR'` : diversification vs expected rate of return
                
            The default is `'RR_risk'`.
        `opt` : optional;
            Additonal parameters:\n
                * `'title'` : `str`; 
                    The default is 'Portfolio frontiers'
                * `'xlabel'` : `str`; 
                    The default is \n
                    - `'risk'` if `fig_type='RR_risk'`,
                    - `'rate of returns'` otherwise.
                * `'ylabel'` : `str`; 
                    The default is \n
                        - `'rate of returns'` if `fig_type='RR_risk'`
                        - `'sharpe'` if `fig_type='Sharpe_RR'`
                        - `'diversification'` if `fig_type=diverse_RR`   
                * `'tangent'` : Boolean; 
                    If set to `True`, then the tangent
                    (to max sharpe point) is added. It has effect only  if
                    `fig_type='RR_risk'`. The default is `True`.
                * `saveto` : `str`; 
                    File name to save the figure. The extension dictates the 
                    format: png, pdf, svg, etc. For more details see the 
                    `mathplotlib` documentation for `savefig`. 
                    The default is `None`.
                * `data` : `defaultdict`; 
                    Numerical data to construct the plot. If it is not `None`,
                    then it will take precedence and no other numerical 
                    evaluations will be performed. 
                    It is meant to produce different plot representations
                    without reevaluations. The default is `None`.
                * `invN_label` : `str`;
                    The label for equal weighted portfolio. 
                    The default is `'1/N'`.
                * `invNrisk_label` : `str`;
                    The lable for efficient portfolio with same risk as 
                    equal weighted portfolio. The default is `invNrisk`.
                * `invNdiverse_label` : `str`;
                    The label for diverse-efficient portfolio with the same
                    diversificeation factor as equal weighted porfolio.
                    The default is `'invNdiv'`.
                * `invNdrr_label` : `str`;
                    The label of diverse-efficient portfolio with the same
                    expected rate of return as equal weighted portolio.
                    The defualt is `'invNdrr'`.
                * `maxdiverse_label` : `str`;
                    The label of maximum diverified portfolio. 
                    The default is `'MaxD'`.
                * `minrisk_label` : `'str'`;
                    The label of minimum risk portfolio.
                    The default is `'MinR'`.
                * `sharpe_label` : `str`;
                    The label of maximum Sharpe portfolio. 
                    The default is `'Sharpe'`.

        Returns
        -------
        `dict`
            Numerical data used to make the plots. It can be passed back to
            reconstruct the plots without reevaluations.
        """
        if fig_type == 'RR_risk':
            xlabel = 'risk'
            ylabel = 'rate of return'
            plot_ff = self._plot_f1
        elif fig_type == 'Sharpe_RR':
            xlabel = 'rate of return'
            ylabel = 'sharpe'
            plot_ff = self._plot_f2
        elif fig_type == 'Diverse_RR':
            xlabel = 'rate of return'
            ylabel = 'diversification'
            plot_ff = self._plot_f3
        else:
            fty = ['RR_risk', 'Sharpe_RR', 'Diverse_RR']
            raise ValueError(f"fig_type={fig_type} ; must be one of {fty}")
                
        options_default = {'invN_label': '1/N',
                           'invNrisk_label': 'InvNrisk',
                           'invNdiverse_label' : 'InvNdiv',
                           'invNdrr_label' : 'InvNdrr',
                           'maxdiverse_label' : 'MaxD',
                           'minrisk_label' : 'MinR',
                           'sharpe_label' : 'Sharpe',
                           'comp_label': True,
                           'addport_label': True,
                           'tangent' : True,
                           'title' : 'Portfolio frontiers',
                           'xlabel' : xlabel,
                           'ylabel' : ylabel,
                           }
        options = defaultdict(lambda: None, options_default)
        for k, v in opt.items():
            if isinstance(v, dict) & (k != 'data'):
                options.update(v)
            else:
                options[k] = v

        resout = options['data']
        del options['data']
 
        # calc
        if resout is None:
            res = defaultdict(lambda: defaultdict(lambda: None))
            
            minrisk_calc = False
            if minrisk or efficient != 0 or inefficient != 0:
                minrisk_calc = True
            
            maxdiverse_calc = False
            self.getRiskComp()
            if maxdiverse or diverse_efficient != 0 or diverse_inefficient != 0: 
                maxdiverse_calc = True
  
            # min risk
            res['risk_min']['risk_min'] = minrisk
            if minrisk_calc:
                ww_min = self.getWeights(rtype='MinRisk')
                if self.status != 0:
                    res['risk_min']['risk_min'] = False
                    warnings.warn("Warning: MinRisk failed")
                else:
                    risk_min = self.risk
                    rr_min = self.RR
                    res['risk_min']['risk'] = risk_min
                    res['risk_min']['rr'] = rr_min
                    res['risk_min']['ww'] = ww_min
                    res['risk_min']['sharpe'] = \
                        (np.dot(ww_min, self.muk) - musharpe) / risk_min \
                            if risk_min != 0 else 0
                    res['risk_min']['diverse'] = self.getDiversification(ww_min)
                    res['risk_min']['label'] = options['minrisk_label']
    
            # efficient frontier
            res['efficient']['efficient'] = 0
            if efficient > 0:
                if rr_min is None:
                    warnings.warn("Warning: MinRisk -> efficient failed")
                    res['efficient']['efficient'] = 0
                else:
                    mu_max = self.muk.max()
                    rr_upper = mu_max - np.abs(mu_max) * _RR_ABS_TOL_
                    rr_eff = np.linspace(rr_min, rr_upper, efficient)
                    risk_eff = []
                    ww_eff = []
                    RR_eff = []
                    sharpe_eff = []
                    diverse_eff = []
                    for rr in rr_eff:
                        _ww = self.getWeights(rtype='Risk', mu=rr)
                        if self.status != 0:
                            warnings.warn(f"Warning: efficient for rr {rr} "
                                          "failed")
                            continue
                        ww_eff.append(_ww)
                        risk_eff.append(self.risk)
                        RR_eff.append(self.RR)
                        sharpe_eff.append(
                            (np.dot(self.ww, self.muk) - musharpe) / self.risk)
                        diverse_eff.append(self.getDiversification(self.ww))
                    res['efficient']['efficient'] = len(risk_eff)
                    res['efficient']['risk'] = np.array(risk_eff)
                    res['efficient']['rr'] = np.array(RR_eff)
                    res['efficient']['ww'] = ww_eff
                    res['efficient']['sharpe'] = np.array(sharpe_eff)
                    res['efficient']['diverse'] = np.array(diverse_eff)
    
            # inefficient frontier
            res['inefficient']['inefficient'] = 0
            if inefficient > 0 :
                if rr_min is None:
                    warnings.warn("Warning: MinRisk -> inefficient failed")
                    res['inefficient']['inefficient'] = 0
                else:
                    mu_min = self.muk.min()
                    rr_lower = mu_min + np.abs(mu_min) * _RR_ABS_TOL_
                    rr_ineff = np.linspace(rr_lower, rr_min, inefficient)
                    risk_ineff = []
                    ww_ineff = []
                    RR_ineff = []
                    sharpe_ineff = []
                    diverse_ineff = []
                    for rr in rr_ineff:
                        _ww = self.getWeights(rtype='Risk', mu=rr, d=-1)
                        if self.status != 0:
                            warnings.warn(f"Warning: inefficient for rr {rr} "
                                          "failed")
                            continue
                        ww_ineff.append(_ww)
                        risk_ineff.append(self.risk)
                        RR_ineff.append(self.RR)
                        sharpe_ineff.append(
                            (np.dot(self.ww, self.muk) - musharpe) / self.risk)
                        diverse_ineff.append(self.getDiversification(self.ww))
                    res['inefficient']['inefficient'] = len(risk_ineff)
                    res['inefficient']['risk'] = np.array(risk_ineff)
                    res['inefficient']['rr'] = np.array(RR_ineff)
                    res['inefficient']['ww'] = ww_ineff
                    res['inefficient']['sharpe'] = np.array(sharpe_ineff)
                    res['inefficient']['diverse'] = np.array(diverse_ineff)
    
            # sharpe point
            res['sharpe']['sharpe'] = sharpe
            if sharpe:
                ww_sharpe = self.getWeights(rtype='Sharpe', mu0=musharpe)
                if self.status != 0:
                    res['sharpe']['sharpe'] = False
                    warnings.warn("Warning: MinRisk failed")
                else:
                    rr_sharpe = self.RR
                    risk_sharpe = self.risk
                    sharpe_sharpe = self.sharpe
                    diverse_sharpe = self.getDiversification(ww_sharpe)
                    res['sharpe']['mu'] = musharpe
                    res['sharpe']['risk'] = risk_sharpe
                    res['sharpe']['rr'] = rr_sharpe
                    res['sharpe']['ww'] = ww_sharpe
                    res['sharpe']['sharpe'] = sharpe_sharpe
                    res['sharpe']['diverse'] = diverse_sharpe
                    res['sharpe']['label'] = options['sharpe_label']
    
            # portfolio components
            res['component']['component'] = component
            if component:
                risk_comp = []
                rr_comp = []
                label_comp = []
                ww_comp = []
                sharpe_comp = []
                for k in range(self.mm):
                    ww = [0.] * self.mm
                    ww[k] = 1.
                    risk_comp.append(self.getRisk(ww))
                    if self.status != 0:
                        warnings.warn("Warning: component "
                                      f"{self.rrate.columns[k]} failed")
                        continue
                    rr_comp.append(self.muk[k])
                    label_comp.append(self.muk.index[k])
                    ww_comp.append(ww)
                    sharpe_comp.append(
                        (self.muk[k] - musharpe) / self.risk \
                            if self.risk != 0 else 0)
                if len(risk_comp) == 0:
                    res['component']['component'] = False
                else:
                    res['component']['risk'] = np.array(risk_comp)
                    res['component']['rr'] = np.array(rr_comp)
                    res['component']['ww'] = ww_comp
                    res['component']['label'] = label_comp \
                        if options['comp_label'] else [None] * self.mm
                    res['component']['sharpe'] = np.array(sharpe_comp)
                    res['component']['diverse'] = np.zeros(self.mm)
    
            # random portfolios
            res['randomport']['randomport'] = 0
            if randomport > 0:
                risk_rp = []
                rr_rp = []
                ww_rp = []
                sharpe_rp = []
                diverse_rp = []
                for _ in range(randomport):
                    ww = self._ww_gen()
                    _div = self.getDiversification(ww)
                    if self.status != 0:
                        warnings.warn(f"Warning: randomport {ww} failed")
                        continue
                    diverse_rp.append(_div)
                    risk_rp.append(self.risk)
                    rr_rp.append(self.RR)
                    ww_rp.append(ww)
                    sharpe_rp.append(
                        (np.dot(self.ww, self.muk) - musharpe) / self.risk)
                res['randomport']['randomport'] = len(risk_rp)
                res['randomport']['risk'] = np.array(risk_rp)
                res['randomport']['rr'] = np.array(rr_rp)
                res['randomport']['ww'] = ww_rp
                res['randomport']['sharpe'] = np.array(sharpe_rp)
                res['randomport']['diverse'] = np.array(diverse_rp)
    
            # invN
            res['invN']['invN'] = invN
            if invN:
                res['invN']['ww'] = np.full(self.mm, 1 / self.mm)
                _div = self.getDiversification(res['invN']['ww'])
                if self.status != 0:
                    res['invN']['invN'] = False
                    warnings.warn("Warning: invN failed")
                else:
                    res['invN']['diverse'] = _div
                    res['invN']['risk'] = self.risk
                    res['invN']['rr'] = self.RR
                    res['invN']['label'] = options['invN_label']
                    res['invN']['sharpe'] = \
                        (np.dot(self.ww, self.muk) - musharpe) / self.risk
                
            # invNrisk
            res['invNrisk']['invNrisk'] = invNrisk
            if invNrisk:
                _ww = self.getWeights(rtype="InvNrisk")
                if self.status != 0:
                    res['invNrisk']['invNrisk'] = False
                    warnings.warn("Warning: invNrisk failed")
                else:
                    res['invNrisk']['ww'] = _ww
                    res['invNrisk']['risk'] = self.risk
                    res['invNrisk']['rr'] = self.RR
                    res['invNrisk']['label'] = options['invNrisk_label']
                    res['invNrisk']['sharpe'] = \
                        (np.dot(self.ww, self.muk) - musharpe) / self.risk
                    res['invNrisk']['diverse'] = \
                        self.getDiversification(res['invNrisk']['ww'])
                    
            # invNdiverse
            res['invNdiverse']['invNdiverse'] = invNdiverse
            if invNdiverse:
                _ww = self.getWeights(rtype="InvNdiverse")
                if self.status != 0:
                    res['invNdiverse']['invNdiverse'] = False
                    warnings.warn("Warning: invNdiverse failed")
                else:
                    res['invNdiverse']['ww'] = _ww
                    res['invNdiverse']['risk'] = self.risk
                    res['invNdiverse']['rr'] = self.RR
                    res['invNdiverse']['label'] = options['invNdiverse_label']
                    res['invNdiverse']['sharpe'] = \
                        (np.dot(self.ww, self.muk) - musharpe) / self.risk
                    res['invNdiverse']['diverse'] = self.diverse
                
            # invNdrr
            res['invNdrr']['invNdrr'] = invNdrr
            if invNdrr:
                _ww = self.getWeights(rtype="InvNdrr")
                if self.status != 0:
                    res['invNdrr']['invNdrr'] = False
                    warnings.warn("Warning: maxdiverse failed")
                else:
                    res['invNdrr']['ww'] = _ww
                    res['invNdrr']['risk'] = self.risk
                    res['invNdrr']['rr'] = self.RR
                    res['invNdrr']['label'] = options['invNdrr_label']
                    res['invNdrr']['sharpe'] = \
                        (np.dot(self.ww, self.muk) - musharpe) / self.risk
                    res['invNdrr']['diverse'] = self.diverse
                
            # maxdiverse
            res['maxdiverse']['maxdiverse'] = maxdiverse
            if maxdiverse_calc:
                ww_d = self.getWeights('MaxDiverse')  
                if self.status != 0:
                    res['maxdiverse']['maxdiverse'] = False
                    warnings.warn("Warning: maxdiverse failed")
                else:
                    res['maxdiverse']['risk'] = self.risk 
                    rr_maxdiverse = self.RR 
                    res['maxdiverse']['rr'] = rr_maxdiverse
                    res['maxdiverse']['ww'] = ww_d
                    res['maxdiverse']['diverse'] = self.diverse
                    res['maxdiverse']['label'] = options['maxdiverse_label']
                    res['maxdiverse']['sharpe'] = \
                        (np.dot(ww_d, self.muk) - musharpe) / self.risk
                
            # diverse-efficient frontier
            res['diverse_efficient']['diverse_efficient'] = 0
            if diverse_efficient > 0 :
                if rr_maxdiverse is None:
                    warnings.warn("Warning: maxdiverse -> diverse_efficient "
                                  "failed")
                    res['diverse_efficient']['diverse_efficient'] = 0
                else:
                    mu_max = self.muk.max() 
                    rr_upper = mu_max - np.abs(mu_max) * _RR_ABS_TOL_
                    rr_d_eff = np.linspace(rr_maxdiverse, rr_upper,
                                           diverse_efficient)
                    risk_d_eff = []
                    ww_d_eff = []
                    RR_d_eff = []
                    sharpe_d_eff = []
                    diverse_d_eff = []
                    for rr in rr_d_eff:
                        _ww = self.getWeights(rtype='Diverse', mu=rr)
                        if self.status != 0:
                            warnings.warn("Warning: diverse_efficient "
                                          f"for rr {rr} failed")
                            continue
                        ww_d_eff.append(_ww)
                        risk_d_eff.append(self.risk)
                        RR_d_eff.append(self.RR)
                        sharpe_d_eff.append(
                            (np.dot(self.ww, self.muk) - musharpe) / self.risk \
                                if self.risk != 0 else 0)
                        diverse_d_eff.append(self.getDiversification(self.ww))
                    res['diverse_efficient']['diverse_efficient'] = len(risk_d_eff)
                    res['diverse_efficient']['risk'] = np.array(risk_d_eff)
                    res['diverse_efficient']['rr'] = np.array(RR_d_eff)
                    res['diverse_efficient']['ww'] = ww_d_eff
                    res['diverse_efficient']['sharpe'] = np.array(sharpe_d_eff)
                    res['diverse_efficient']['diverse'] = np.array(diverse_d_eff)
                
            # diverse-inefficient frontier
            res['diverse_inefficient']['diverse_inefficient'] = 0
            if diverse_inefficient > 0 :
                if rr_maxdiverse is None:
                    warnings.warn("Warning: maxdiverse -> diverse_inefficient "
                                  "failed")
                    res['diverse_inefficient']['diverse_inefficient'] = 0
                else:
                    rr_min = self.muk.min()
                    rr_lower = rr_min + np.abs(mu_min) * _RR_ABS_TOL_
                    rr_d_ineff = np.linspace(rr_lower, rr_maxdiverse, 
                                             diverse_inefficient)
                    risk_d_ineff = []
                    ww_d_ineff = []
                    RR_d_ineff = []
                    sharpe_d_ineff = []
                    diverse_d_ineff = []
                    for rr in rr_d_ineff:
                        _ww = self.getWeights(rtype='Diverse', mu=rr, d=-1)
                        if self.status != 0:
                            warnings.warn(f"Warning: diverse_inefficient "
                                          f"for rr {rr} failed")
                            continue
                        ww_d_ineff.append(_ww)
                        risk_d_ineff.append(self.risk)
                        RR_d_ineff.append(self.RR)
                        sharpe_d_ineff.append(
                            (np.dot(self.ww, self.muk) - musharpe) / self.risk \
                                if self.risk != 0 else 0)
                        diverse_d_ineff.append(self.getDiversification(self.ww))
                    res['diverse_inefficient']['diverse_inefficient'] = \
                        len(risk_d_ineff)
                    res['diverse_inefficient']['risk'] = np.array(risk_d_ineff)
                    res['diverse_inefficient']['rr'] = np.array(RR_d_ineff)
                    res['diverse_inefficient']['ww'] = ww_d_ineff
                    res['diverse_inefficient']['sharpe'] = \
                        np.array(sharpe_d_ineff)
                    res['diverse_inefficient']['diverse'] = \
                        np.array(diverse_d_ineff)
                
            # addport
            res['addport']['addport'] = False
            if addport is not None:
                if isinstance(addport, dict):
                    aport = pd.DataFrame().from_dict(addport, 'index', 
                                                    columns=self.rrate.columns)
                elif isinstance(addport, pd.core.frame.DataFrame):
                    aport = addport[self.rrate.columns]
                else:
                    raise ValueError("addport must be a dict or a "
                                     + "pandas.DataFrame, see doc")
     
                def ff_addport(x):
                    diverse_addport = self.getDiversification(x)
                    if self.status != 0:
                        warnings.warn(f"Warning: addport {x} failed")
                        risk_addport = None
                        RR_addport = None
                        sharpe_addport = None
                    else:
                        risk_addport = self.risk
                        RR_addport = self.RR
                        sharpe_addport = (np.dot(x, self.muk) - musharpe) \
                             / risk_addport if risk_addport != 0  else 0 
                             
                    return pd.Series([self.status, risk_addport, RR_addport, 
                        sharpe_addport, diverse_addport],
                        index=['status','risk', 'rr', 'sharpe', 'diverse'])
                
                res_addport = aport.apply(lambda x: ff_addport(x), axis=1).dropna()
                if res_addport.size > 0:
                    res['addport']['addport'] = True
                    res['addport']['risk'] = res_addport['risk'].to_list()
                    res['addport']['rr'] = res_addport['rr'].to_list()
                    res['addport']['sharpe'] = res_addport['sharpe'].to_list()
                    res['addport']['diverse'] = res_addport['diverse'].to_list()
                    res['addport']['ww'] = aport.to_numpy().tolist()
                    res['addport']['label'] = \
                        aport.index.to_list() if options['addport_label'] \
                            else [None] * aport.shape[0]
                    res['addport']['hidden_label'] = aport.index.to_list()
                
            resout = res
            res['options'] = options
        else:
            res = copy.deepcopy(resout)

            res['options'] = options
            res['risk_min']['label'] = options['minrisk_label']
            res['sharpe']['label'] = options['sharpe_label']
            res['component']['label'] = self.muk.index.tolist() \
                if options['comp_label'] else [None] * self.mm
            res['invN']['label'] = options['invN_label']
            res['invNrisk']['label'] = options['invNrisk_label']
            res['invNdiverse']['label'] = options['invNdiverse_label']
            res['invNdrr']['label'] = options['invNdrr_label']
            res['maxdiverse']['label'] = options['maxdiverse_label']
            if res['addport']['addport']:
                if options['addport_label']:
                    res['addport']['label'] = res['addport']['hidden_label']
                else:
                    res['addport']['label'] = [None] * \
                        len(res['addport']['hidden_label'])
            
        # Plot
        plot_ff(res)

        return resout


    def _plot_f1(self, res):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.set(title=res['options']['title'], 
               xlabel=res['options']['xlabel'], 
               ylabel=res['options']['ylabel'])

        if res['risk_min']['risk_min']:
            ax.scatter(x=res['risk_min']['risk'], y=res['risk_min']['rr'],
                       marker='D', c='g', s=36)
            ax.annotate(res['risk_min']['label'],
                        (res['risk_min']['risk'] * 1.01,
                         res['risk_min']['rr'] * 1.04))

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

        for invstr in ['invN', 'invNrisk', 'invNdiverse', 'invNdrr']:
            if res[invstr][invstr]:
                ax.scatter(x=res[invstr]['risk'], y=res[invstr]['rr'],
                          marker='x', c='g', s=36)
                lf = 1.01
                ax.annotate(res[invstr]['label'],
                            (res[invstr]['risk'] * lf,
                             res[invstr]['rr'] * lf))
  
        if res['maxdiverse']['maxdiverse']:
            ax.scatter(x=res['maxdiverse']['risk'], y=res['maxdiverse']['rr'],
                       marker='x', c='b', s=36)
            lf = 1.01
            ax.annotate(res['maxdiverse']['label'],
                        (res['maxdiverse']['risk'] * lf,
                         res['maxdiverse']['rr'] * lf))
    
        if res['diverse_efficient']['diverse_efficient'] > 0:
            ax.plot(res['diverse_efficient']['risk'], 
                    res['diverse_efficient']['rr'],
                    color='lightblue', linewidth=1, linestyle=':')

        if res['diverse_inefficient']['diverse_inefficient'] > 0:
            ax.plot(res['diverse_inefficient']['risk'], 
                    res['diverse_inefficient']['rr'],
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
            ax.annotate(res['sharpe']['label'],
                        (res['sharpe']['risk'] * 0.91,
                         res['sharpe']['rr'] * 1.04))
            if res['options']['tangent']:
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
                    
        if res['options']['saveto'] is not None:
            plt.savefig(res['options']['saveto'])
        plt.show()


    def _plot_f2(self, res):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.set(title=res['options']['title'], 
               xlabel=res['options']['xlabel'], 
               ylabel=res['options']['ylabel'])

        if res['risk_min']['risk_min']:
            ax.scatter(x=res['risk_min']['rr'], y=res['risk_min']['sharpe'], 
                       marker='D', c='g', s=36)
            ax.annotate(res['risk_min']['label'],
                        (res['risk_min']['rr'] * 0.92,
                         res['risk_min']['sharpe'] * 1.04))

        if res['efficient']['efficient'] > 0:
            ax.plot(res['efficient']['rr'], res['efficient']['sharpe'],
                    color='lightblue', linewidth=2)

        if res['inefficient']['inefficient'] > 0:
            ax.plot(res['inefficient']['rr'], res['inefficient']['sharpe'],
                    color='red', linewidth=2)

        if res['component']['component']:
            lf = 1.01
            for k in range(len(res['component']['risk'])):
                ax.scatter(x=res['component']['rr'][k], 
                           y=res['component']['sharpe'][k],
                           marker='s', c='b', s=16)
                ax.annotate(res['component']['label'][k],
                            (res['component']['rr'][k] * lf, 
                             res['component']['sharpe'][k] * lf))
            
        if res['sharpe']['sharpe']:
            ax.scatter(x=res['sharpe']['rr'], y=res['sharpe']['sharpe'],
                       marker='D', c='g', s=26)
            ax.annotate(res['sharpe']['label'],
                        (res['sharpe']['rr'] * 0.93,
                         res['sharpe']['sharpe'] * 1.02))

        if res['randomport']['randomport'] > 0:
            ax.scatter(x=res['randomport']['rr'], 
                       y=res['randomport']['sharpe'],
                       marker='o', s=9, edgecolors='b', alpha=0.5 )

        for invstr in ['invN', 'invNrisk', 'invNdiverse', 'invNdrr']:
            if res[invstr][invstr]:
                ax.scatter(x=res[invstr]['rr'], y=res[invstr]['sharpe'],
                          marker='x', c='g', s=36)
                lf = 1.01
                ax.annotate(res[invstr]['label'],
                            (res[invstr]['rr'] * lf,
                             res[invstr]['sharpe'] * lf))
   
        if res['maxdiverse']['maxdiverse']:
            ax.scatter(x=res['maxdiverse']['rr'], 
                       y=res['maxdiverse']['sharpe'],
                       marker='x', c='b', s=36)
            lf = 1.01
            ax.annotate(res['maxdiverse']['label'],
                       (res['maxdiverse']['rr'] * lf, 
                        res['maxdiverse']['sharpe'] * lf))

        if res['diverse_efficient']['diverse_efficient'] > 0:
            ax.plot(res['diverse_efficient']['rr'], 
                    res['diverse_efficient']['sharpe'],
                    color='lightblue', linewidth=1, linestyle=':')

        if res['diverse_inefficient']['diverse_inefficient'] > 0:
            ax.plot(res['diverse_inefficient']['rr'], 
                    res['diverse_inefficient']['sharpe'],
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

        if res['options']['saveto'] is not None:
            plt.savefig(res['options']['saveto'])
        plt.show()


    def _plot_f3(self, res):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set(title=res['options']['title'], 
               xlabel=res['options']['xlabel'], 
               ylabel=res['options']['ylabel'])

        if res['risk_min']['risk_min']:
            ax.scatter(x=res['risk_min']['rr'], y=res['risk_min']['diverse'], 
                       marker='D', c='g', s=36)
            lf = 0.95
            ax.annotate(res['risk_min']['label'],
                       (res['risk_min']['rr'] * lf, 
                        res['risk_min']['diverse'] * lf))
    
        if res['efficient']['efficient'] > 0:
            ax.plot(res['efficient']['rr'], res['efficient']['diverse'],
                    color='lightblue', linewidth=1, linestyle=':')
    
        if res['inefficient']['inefficient'] > 0:
            ax.plot(res['inefficient']['rr'], res['inefficient']['diverse'],
                    color='red', linewidth=1, linestyle=':')
    
        if res['component']['component']:
            lf = 1.01
            for k in range(len(res['component']['risk'])):
                ax.scatter(x=res['component']['rr'][k], 
                           y=res['component']['diverse'][k],
                           marker='s', c='b', s=16)
                ax.annotate(res['component']['label'][k],
                            (res['component']['rr'][k] * lf, 
                             res['component']['diverse'][k] * lf))
    
        if res['sharpe']['sharpe']:
            ax.scatter(x=res['sharpe']['rr'], 
                       y=res['sharpe']['diverse'],
                       marker='D', c='g', s=26)
            lf = 1.01
            ax.annotate(res['sharpe']['label'],
                       (res['sharpe']['rr'] * lf, 
                        res['sharpe']['diverse'] * lf))
    
        if res['randomport']['randomport'] > 0:
            ax.scatter(x=res['randomport']['rr'], 
                       y=res['randomport']['diverse'],
                       marker='o', s=9, edgecolors='b', alpha=0.5 )
    
        for invstr in ['invN', 'invNrisk', 'invNdiverse', 'invNdrr']:
            if res[invstr][invstr]:
                ax.scatter(x=res[invstr]['rr'], 
                           y=res[invstr]['diverse'],
                           marker='x', c='g', s=36)
                lf = 1.01
                ax.annotate(res[invstr]['label'],
                            (res[invstr]['rr'] * lf,
                             res[invstr]['diverse'] * lf))
 
        if res['maxdiverse']['maxdiverse']:
            ax.scatter(x=res['maxdiverse']['rr'], 
                       y=res['maxdiverse']['diverse'],
                       marker='x', c='b', s=36)
            lf = 1.01
            ax.annotate(res['maxdiverse']['label'],
                       (res['maxdiverse']['rr'] * lf, 
                        res['maxdiverse']['diverse'] * lf))
            
        if res['diverse_efficient']['diverse_efficient'] > 0:
            ax.plot(res['diverse_efficient']['rr'], 
                    res['diverse_efficient']['diverse'],
                    color='lightblue', linewidth=2)
    
        if res['diverse_inefficient']['diverse_inefficient'] > 0:
            ax.plot(res['diverse_inefficient']['rr'], 
                    res['diverse_inefficient']['diverse'],
                    color='red', linewidth=2)
            
        if res['addport']['addport']:
            lf = 1.01
            for k in range(len(res['addport']['risk'])):
                ax.scatter(x=res['addport']['rr'][k],
                           y=res['addport']['diverse'][k],
                           marker='o', s=9, edgecolors='r', alpha=0.5)
                ax.annotate(res['addport']['label'][k],
                            (res['addport']['rr'][k] * lf,
                             res['addport']['diverse'][k] * lf))
            
        xl1, xl2 = plt.xlim()
        yl1, yl2 = plt.ylim()
        lff = 1.02
        xl2 *= lff
        yl2 *= lff
        ax.set(xlim=(xl1, xl2), ylim=(yl1, yl2))
    
        if res['options']['saveto'] is not None:
            plt.savefig(res['options']['saveto'])
        plt.show()


    def set_random_seed(self, seed = 42):
        """
        Sets the seed for Dirichlet random generator used in viewFrontiers
        to select the random inefficient portfolios.

        Parameters
        ----------
        `seed` : `int`, optional;
            The random generator seed, in case you want to set it to a weird
            value other than `42` :). The default is `42`.

        Returns
        -------
        `None`
        """
        self.rng = np.random.RandomState(seed)


    def _ww_gen(self):
        return self.rng.dirichlet([0.5] * self.mm)
    
    
    def _set_lp_method(self, method):
        self.methods = ['ecos', 'highs-ds', 'highs-ipm', 'highs',
                       'interior-point', 'glpk', 'cvxopt']
        if not method in self.methods:
            raise ValueError(f"unknown method {method} - "
                             f"must be one of {self.methods}")
        self.method = method
        
        
    def _set_qp_method(self, method):
        self.methods = ['ecos', 'cvxopt']
        if not method in self.methods:
            raise ValueError(f"unknown method {method} - "
                             f"must be one of {self.methods}")
        self.method = method
        
        
    def _set_socp_method(self, method):
        self.methods = ['ecos', 'cvxopt']
        if not method in self.methods:
            raise ValueError(f"unknown method {method} - "
                             f"must be one of {self.methods}")
        self.method = method
        

    
    def getRiskComp(self):
        """
        Returns the risk of each portfolio component.

        Returns
        -------
        `pandas.Series`
        """
        if self._flag_risk_comp_calc: 
            self.status = self._status_risk_comp_calc
            return self.risk_comp
            
        srisk = np.zeros(self.rrate.shape[1])
        status = 0
        for i in range(len(srisk)):
            wws = np.zeros(len(srisk))
            wws[i] = 1
            srisk[i] = self.getRisk(wws)
            if self.status != 0:
                status = self.status
                srisk[i] = None
 
        self.risk_comp = pd.Series(srisk, self.rrate.columns)
        self._flag_risk_comp_calc = True
        self.status = status
        self._status_risk_comp_calc = self.status
        
        if self.status != 0:
            warnings.warn(f"Fail to compute for edate {self.rrate.index[-1]} "
                          f"all risk coml\n{self.risk_comp}")

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
            be in order of `rrate.columns`. If it is a `pandas.Series` then 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (not necessary in the same order).
        `rrate` : `pandas.DataFrame`, optional;
            Contains the portfolio components historical
            rates of returns. If it is not `None` then it will overwrite the
            `rrate` computed in the constructor from `mktdata`.
            The default is `None`.

        Returns
        -------
        `float` : The diversification value. \n
        """
        self._reset_output()
        
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
        if self.status != 0:
            warnings.warn("Warning: Fail to compute all risk comp\n"
                          f"{self.risk_comp}")
            return None
            
        wrisk = self.getRisk(w)
        if self.status !=0:
            warnings.warn(f"Warning: Fail to compute risk for {ww}")
            return None
        
        crisk = np.dot(self.risk_comp, w)

        if crisk == 0:
            self.diverse = 0
        else:
            self.diverse = 1 - wrisk / crisk
        
        return self.diverse
    
    
    def set_method(self, method):
        """
        Set computation method

        Parameters
        ----------
        `method` : `str`;
            Must be a valid method name. It will overright the value 
            set by the constructor.

        Returns
        -------
        `None`
        """
        self._set_method(method)
    
    
    def _norm_weights(self):
        self.ww[self.ww < _WW_TOL_] = 0.
        self.ww /= self.ww.sum()
        
        
    def _single_asset_port(self, port):
        if port == 'min':
            psymb = self.muk.idxmin()
        elif port == 'max':
            psymb = self.muk.idxmax()
        elif port in self.rrate.columns:
            psymb = port
        else:
            raise ValueError(f"argumkent ({port}) must be "
                             f"one of 'min', 'max', {self.symb}")
        w = pd.Series(0., index=self.rrate.columns)
        w[psymb] = 1.
        self.getRisk(w) 
        self.diverse = 0.


    # to be implemented in the deriv classes
    def _risk_calc(self, prate, alpha):
        # Computes a simple risk measure for `prate` time series and `alpha`
        # parameter (where is the case)
        pass
    def _risk_min(self, d=1):
        # Computes the optimal portfolio with minimum risk for fixed expected 
        # rate of return where `d` define the risk frontier type: 
        # +1 for efficient and -1 for inefficient 
        pass
    def _sharpe_max(self):
        # Computes the optimal portfolio with maximum Sharpe ratio
        pass
    def _sharpe_inv_min(self):
        # Computes the optimal portfolio with minimum inverse Sharpe ratio
        pass
    def _rr_max(self):
        # Computes the optimal-risk portfolio with maximum expected rate of 
        # return for fixed risk
        pass
    def _risk_averse(self):
        # Computes the optimal-risk portfolio for fixed risk-aversion factor
        pass
    def _risk_diversification(self, d=1):
        # Computes the optimal-diversified portfolio for fixed expected
        # rate of return where `d` define the diversified frontier type: 
        # +1 for efficient and -1 for inefficient.
        # Minimum of 1-D.
        pass
    def _risk_inv_diversification(self, d=1):
        # Computes the optimal diversified portfolio for fixed expected
        # rate of return where `d` define the diversified frontier type: 
        # +1 for efficient and -1 for inefficient.
        # Maximum of inverse 1-D. 
        pass
    def _rr_max_diversification(self):
        # Computes the optimal diversified portfolio for fixed diversification
        # factor. Maximization of expected rate of return.
        pass
    def _set_method(self, method):
        # sets method
        pass
    