import numpy as np
import scipy.sparse as sps
import time

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _lp_solver

class GINIAnalyzer(_RiskAnalyzer):
    """
    GINI based optimal portfolio strategies.
    
    **Attributes**
        * `status` : `int` - the computation status (`0` - success, 
          any other value signifies an error)
        * `ww` : `pandas.Series` -  the portfolio weights 
        * `RR` : `float` - portfolio rate of return
        * `risk` : `float` - portfolio GINI risk
        * `primary_risk_comp` : `list` - redundant (single element list 
          containing GINI risk value)
        * `secondary_risk_comp` : `list` - redundant 
          (same as `primary_risk_comp`)
        * `sharpe` : `float` - GINI-Sharpe ration if `rtype` is set to 
          `'Sharpe'` or `'Sharpe2'` otherwise `None`. 
        * `diverse` : `float` - diversification factor if `rtype` is set 
          to `'Divers'` or `'MaxDivers'` otherwise `None`.
        * `name` : `str` - portfolio name
        
    Note the following 2 important methods:
        * `getWeights` : Computes the optimal portfolio weights.
          During its computations the following class members are also set:
          `risk`, `primary_risk_comp`, `secondary_risk_comp`, `sharpe`,  `RR`, 
          `divers`.
        * `getPositions` : Provides practical information regarding the portfolio
          rebalancing delta positions and costs.  
    """
    def __init__(self, mktdata=None, colname='adjusted', freq='Q', 
                 hlength=1.25, name='GINI', rtype='Sharpe', mu=None, 
                 d=1, mu0=0., aversion=None, ww0=None, method='ecos',
                 verbose=False):
        """
        Constructor

        Parameters
        ----------
        mktdata : `pandas.DataFrame`, optional
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        colname : `str`, optional
            Name of the price column from mktdata used in the weight's 
            calibration. The default is `'adjusted'`.
        freq : `str`, optional
            Rate of return horizon. It could be 
            `'Q'` for a quarter or `'M'` for a month. The default is `'Q'`.
        hlength : `float`, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is `1.25` years.
        name : `str`, optional
            Portfolio name. The default is `'GINI'`.
        rtype : `str`, optional
            Optimization type. Possible values: \n
                `'Risk'` : optimal risk portfolio for targeted expected rate of 
                return.\n
                `'Sharpe'` : optimal Sharpe portfolio - maximization solution.\n
                `'Sharpe2'` : optimal Sharpe portfolio - minimization solution.\n
                `'MinRisk'` : minimum risk portfolio.\n
                `'RiskAverse'` : optimal risk portfolio for a fixed 
                risk-aversion factor.\n
                `'InvNrisk'` : optimal risk portfolio with the same risk value 
                as a benchmark portfolio (e.g., same as equal weighted 
                portfolio).\n
                `'Diverse'` : optimal diversified portfolio for targeted
                expected rate of return (max of inverse 1-Diverse).\n
                `'Diverse2'` : optimal diversified portfolio for targeted
                expected rate of return (min of 1-Diverse).\n
                `'MaxDiverse'` : maximum diversified portfolio.\n
                `'InvNdiverse'` : optimal diversified portfolio with the same
                diversification factor as a benchmark portfolio 
                (e.g., same as equal weighted portfolio).\n
                `'InvNdrr'` : optimal diversified portfolio with the same 
                expected rate of return as a benchmark portfolio
                (e.g., same as equal weighted portfolio).\n
            The default is `'Sharpe'`.
        mu : `float`, optional
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype='Risk'` or `rtype='Divers'`.
            The default is `None`.
        d : `int`, optional
            Frontier type. Active only if `rtype='Risk'`. A value of `1` will
            trigger the evaluation of optimal portfolio along the efficient
            frontier. Otherwise, it will find the portfolio with the lowest
            rate of return along the inefficient portfolio frontier.
            The default is `1`.
        mu0 : `float`, optional
            Risk-free rate accessible to the investor.
            Relevant only if `rtype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        aversion : `float`, optional
            The value of the risk-aversion coefficient.
            Must be positive. Relevant only if `rtype='RiskAverse'`.
            The default is `None`.
        ww0 : `list`, `numpy.array` or `pandas.Series`, optional
            Targeted portfolio weights. 
            Relevant only if `rtype='InvNrisk'`.
            Its length must be equal to the number of symbols in `rrate` 
            (mktdata). All weights must be >= 0 with sum > 0.
            If it is a `list` or a `numpy.array` then the weights are assumed 
            to be in order of `rrate.columns`. If it is a `pandas.Series` then 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (same symbols, not necessarily in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.
        method : `str`, optional
            Linear programming numerical method. 
            Could be: `'ecos'`, `'highs-ds'`, `'highs-ipm'`, `'highs'`, 
            `'interior-point'`, `'glpk'` and `'cvxopt'`.
            The default is `'ecos'`.
        verbose : `Boolean`, optional
            If it is set to `True`, then various computation messages 
            (meant as warnings) will be printed. The default is `False`.
          
        Returns
        -------
        The object.
        """
        self.drate = None
        self.nn2 = None
        super().__init__(mktdata, colname, freq, hlength, name,
                         rtype, mu, d, mu0, aversion, ww0, verbose)
        self._set_method(method)
        
        
    def _set_method(self, method):
        self._set_lp_method(method)
        
    
    def set_rrate(self, rrate):
        """
        Sets the MkT Data.

        Parameters
        ----------
        rrate : `pandas.DataFrame`
            Market data. It will overwrite the value set by the constructor.
        Returns
        -------
        `None`
        """
        if rrate is None:
            # nothing to do
            return
        
        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean(numeric_only=True)
        self.rrate = rrate.copy()
        
        self.nn2 = int(self.nn * (self.nn - 1) / 2)
        yy = []
        for m in range(self.mm):
            x = self.rrate.iloc[:,m]
            yy.append([x.iloc[i] - x.iloc[j] \
                       for i in range(self.nn - 1) \
                       for j in range(i + 1, self.nn)])
        self.drate = np.concatenate(yy)
        
        
    def _risk_calc(self, prate, alpha):
        nn = len(prate)
 
        y = np.sort(prate, kind='heapsort')
        gini = np.sum([y[i] * ( 2 * i - nn + 1) for i in range(nn)]) / nn / nn
        # status, gini, gini
        return 0, gini, gini
    
    
    def _risk_min(self, d=1):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # in total dim = mm + nn2 
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = [0.] * mm + [norm] * nn2
        
        # bild G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += [nn2 * 2] * mm
        G_data += list(-self.muk * d)
        
        G_icol += list(range(mm + nn2))
        G_irow += list(range(nn2 * 2 + 1, nn2 * 3 + 1 + mm))
        G_data += [-1.] * (mm + nn2)
        
        G_shape = (nn2 * 3 + 1 + mm, mm + nn2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
   
        # build h
        h_data = [0.] * (nn2 * 2) + [-self.mu * d] + [0.] * (mm + nn2)
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + nn2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # GINI
        self.risk = res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # defalut to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _sharpe_max(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # t <- mm + nn2
        # in total dim = mm + nn2 + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = list(-self.muk) + [0.] * nn2 + [self.mu]
        
        # bild G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * mm
        
        G_shape = (nn2 * 2 + mm , mm + nn2 + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm) 
        
        # build A
        A_icol = list(range(mm, mm + nn2)) + list(range(mm)) + [mm + nn2] 
        A_irow = [0] * nn2 + [1] * (mm + 1)
        A_data = [norm] * nn2 + [1.] * mm + [-1.]
        
        A_shape = (2, mm + nn2 + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Gini-Sharpe
        self.sharpe = -res['pcost']
        # GINI
        self.risk = 1. / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost'] / t + self.mu
        # default to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _sharpe_inv_min(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # t <- mm + nn2
        # in total dim = mm + nn2 + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = [0.] * mm + [norm] * nn2 + [0.]
        
        # build G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * mm
        
        G_shape = (nn2 * 2 + mm, mm + nn2 + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm) 
        
        # build A
        A_icol = (list(range(mm)) + [mm + nn2]) * 2
        A_irow = [0] * (mm + 1) + [1] * (mm + 1)
        A_data = list(self.muk) + [-self.mu] + [1.] * mm + [-1.]
        
        A_shape = (2, mm + nn2 + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Gini-Sharpe
        self.sharpe = 1. / res['pcost']
        # GINI
        self.risk = res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # default to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _rr_max(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # in total dim = mm + nn2
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = list(-self.muk) + [0.] * nn2
        
        # build G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * mm
        
        G_shape = (nn2 * 2 + mm, mm + nn2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm) 
        
        # build A
        A_icol = list(range(mm)) + list(range(mm, mm + nn2))
        A_irow = [0] * mm + [1] * nn2
        A_data = [1.] * mm + [norm] * nn2
        
        A_shape = (2, mm + nn2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., self.risk]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost']
        # default to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww


    def _risk_averse(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # in total dim = mm + nn2 
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = list(-self.muk) + [self.Lambda * norm] * nn2
        
        # bild G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * (mm)
        
        G_shape = (nn2 * 2 + mm, mm + nn2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm)
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_shape = (1, mm + nn2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # GINI
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # default to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _risk_diversification(self, d=1):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # t <- mm + nn2
        # in total dim = mm + nn2 + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = [0.] * mm + [norm] * nn2 + [0.]
        
        # build G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm)) + [mm + nn2]
        G_irow += [nn2 * 2] * (mm + 1)
        G_data += list(-self.muk * d) + [self.mu * d]
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2 + 1, nn2 * 2 + mm + 1))
        G_data += [-1.] * mm
        
        G_shape = (nn2 * 2 + mm + 1, mm + nn2 + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm + 1) 
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = list(self.risk_comp)
        
        A_icol += list(range(mm)) + [mm + nn2]
        A_irow += [1] * (mm + 1)
        A_data += [1.] * mm + [-1.]
        
        A_shape = (2, mm + nn2 + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
           if self.verbose:
               print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                     f"status {res['status']} :: {res['infostring']}")
           return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Gini-Divers
        self.diverse= 1. - res['pcost']
        # Gini
        self.risk = res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # default to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _risk_inv_diversification(self, d=1):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # t <- mm + nn2
        # in total dim = mm + nn2 + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = list(-self.risk_comp) + [0.] * (nn2 + 1)
        
        # build G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm)) + [mm + nn2]
        G_irow += [nn2 * 2] * (mm + 1)
        G_data += list(-self.muk * d) + [self.mu * d]
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2 + 1, nn2 * 2 + mm + 1))
        G_data += [-1.] * mm
        
        G_shape = (nn2 * 2 + mm + 1, mm + nn2 + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm + 1) 
        
        # build A
        A_icol = list(range(mm, mm + nn2))
        A_irow = [0] * nn2
        A_data = [norm] * nn2
        
        A_icol += list(range(mm)) + [mm + nn2]
        A_irow += [1] * (mm + 1)
        A_data += [1.] * mm + [-1.]
        
        A_shape = (2, mm + nn2 + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Gini-Divers
        self.diverse= 1. + 1. / res['pcost']
        # Gini
        self.risk = 1. / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # default to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
        

    def _rr_max_diversification(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # in total dim = mm + nn2
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = list(-self.muk) + [0.] * nn2
        
        # build G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * mm
        
        G_shape = (nn2 * 2 + mm, mm + nn2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm) 
        
        # build A
        A_icol = list(range(mm)) + list(range(mm, mm + nn2)) + list(range(mm))
        A_irow = [0] * mm + [1] * (nn2 + mm)
        A_data = [1.] * mm + [norm] * nn2 + list((self.diverse- 1.) * self.risk_comp)
        
        A_shape = (2, mm + nn2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost']
        # risk
        self.risk = norm * np.sum(res['x'][mm:])
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        # diversification
        self.diverse= 1. - self.risk / np.dot(self.ww, self.risk_comp)
        
        return self.ww