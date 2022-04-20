import numpy as np
import scipy.sparse as sps
import warnings

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _lp_solver

class OmegaAnalyzer(_RiskAnalyzer):
    """
    Omega measure/ratio based portfolio optimization.
    
    Methods:
        * getWeights
        * getRisk
        * getPositions
        * viewForntiers
        * set_rrate
        * set_mktdata
        * set_rtype
        * set_random_seed
    """
    def __init__(self, mu0 = 0., 
                 mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, calendar=None, 
                 rtype='Sharpe', method='ecos'):
        """
        Constructor

        Parameters
        ----------
        mu0 : float, optional
            Omega threshold. The default is 0.
        mktdata : pandas.DataFrame, optional
            Historic daily market data for portfolio components in the format
            returned by azapy.mktData function. The default is None.
        colname : str, optional
            Name of the price column from mktdata used in the weights 
            calibration. The default is 'adjusted'.
        freq : str, optional
            Rate of returns horizon in number of business day. it could be 
            'Q' for quarter or 'M' for month. The default is 'Q'.
        hlength : float, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is 3.25 years.
        calendar : numpy.busdaycalendar, optional
            Business days calendar. If is it `None` then the calendar will be set
            to NYSE business calendar.
            The default is `None`.
        rtype : str, optional
            Optimization type. Possible values \n
                "Risk" : minimization of dispersion (risk) measure for a fixed 
                vale of expected rate of return. \n
                "Sharpe" : maximization of generalized Sharpe ratio.\n
                "Sharpe2" : minimization of the inverse generalized Sharpe 
                ratio.\n
                "MinRisk" : optimal portfolio with minimum dispersion (risk) 
                value.\n
                "InvNRisk" : optimal portfolio with the same dispersion (risk)
                value as equal weighted portfolio. \n
                "RiskAverse" : optimal portfolio for a fixed value of risk 
                aversion coefficient.
            The default is "Sharpe". 
        method : string, optional
            Linear programming numerical method. 
            Could be: 'ecos', 'highs-ds', 'highs-ipm', 'highs', 
            'interior-point', 'glpk' and 'cvxopt'.
            The defualt is 'ecos'.

        Returns
        -------
        The object.
        """
        super().__init__(mktdata, colname, freq, hlength, calendar, rtype)
        
        lp_methods = ['ecos', 'highs-ds', 'highs-ipm', 'highs', 
                       'interior-point', 'glpk', 'cvxopt']
        if not method in lp_methods:
            raise ValueError(f"method must be one of {lp_methods}")
        self.method = method
        
        self.alpha = [mu0]
        
    def viewFrontiers(self, efficient=20, inefficient=20, musharpe=None,
                      component=True, randomport=20, fig_type='RR_risk',
                      options=None, saveto=None, data=None):
        """
        Computes the elements of the portfolio frontiers.

        Parameters
        ----------
        efficient : int, optional
            Number of points along the optimal frontier (equally spaced along 
            the rate of returns). The default is 20.
        inefficient : int, optional
            Number of points along the inefficient frontier (equally spaced 
            along the rate of returns). The default is 20.
        musharpe : float, optional
            Value for the risk-free rate of return used in the evaluation of
            generalized Sharpe ratio. The default is 0.
        component : Boolean, optional
            If True the portfolios containing a single component are evaluated 
            and added to the plot for references. The default is True.
        randomport : int, optional
            The number of portfolios with random weights (inefficient) to be 
            evaluate and added to the plot for reference. The default is 20.
        inverseN : Boolean, optional
            If `True` the equally weighted portfolio and the optimal portfolio 
            with the same dispersion (risk) value are evaluated and added to 
            the plot. The default is `True`.
        fig_type : string, optional
            Graphical representation format.
            If it is set to "RR_risk" the data is plotted in the rate of return 
            vs dispersion representation, otherwise the Sharpe vs rate of 
            return will be used. The default is 'RR_risk'.
        options : dict, optional
            Additional graphical setups (keys): \n
            "title", "xlabel", "ylabel" : strings overwriting the default 
            values. \n
            "tangent" : boolean. If set to `True` it will print
            the Sharpe tangent. The default is `True`.
        save : str, optional
            File name to save the plot. The default is None.
        data : dict, optional
            Numerical data to construct the plot. If it is not `None` it 
            will take precedence and no other numerical evaluation will be 
            performed. It is meant to produce different plot representations
            without recomputation. The default is `None`.

        Returns
        -------
        dict
            Numerical data used to make the plots. 
        """
        if musharpe is not None:
            self.alpha[0] = musharpe
            
        return \
        super().viewFrontiers(efficient=efficient, inefficient=inefficient,
                              musharpe=self.alpha[0],
                              component=component, randomport=randomport,
                              fig_type=fig_type, options=options, 
                              saveto=saveto, data=data)
        
    def set_rrate(self, rrate):
        """
        Sets portfolio components historical rates of returns in the format 
        "date", "symbol1", "symbol2", etc. 

        Parameters
        ----------
        rrate : pandas.DataFrame
            Portfolio components historical rates of returns in the format 
           "date", "symbol1", "symbol2", etc.
            It will overwrite the values set by the constructor.
        Returns
        -------
        None.
        """
        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean()
        self.rrate = rrate
        
    def _risk_calc(self, prate, alpha):
        rr = alpha - prate
        rr[rr < 0] = 0.
        rho = np.mean(rr)
        # status, rho, rho
        return 0, rho, rho
    
    def _risk_min(self, d=1):
        # Order of variables
        # w <- [0 : mm]
        # s <- [mm : mm + nn]
        # in total dim = mm + nn
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = [0.] * mm + [1. / nn] * nn
        
        # build G
        G_icol = list(range(mm)) * nn + list(range(mm, mm + nn)) \
            + list(range(mm))
        G_irow = [k  for k in range(nn) for _ in range(mm)] \
            + list(range(nn)) + [nn] * mm
        G_data = list(np.ravel(-self.rrate)) + [-1.] * nn + list(-self.muk * d)
        G_icol += list(range(mm + nn))
        G_irow += list(range(nn + 1, nn + 1 + mm + nn))
        G_data += [-1.] *(mm + nn)
        
        G_shape = (nn + 1 + mm + nn, mm + nn)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
       
        # build h
        h_data = [-self.alpha[0]] * nn + [-self.mu * d] + [0.] *(mm + nn)
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + nn)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
       
        # build b
        b_data = [1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # delta-risk
        self.risk = res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _sharpe_inv_min(self):
        # Order of variables:
        # w <- [0:mm] 
        # s <- [mm : mm + nn]
        # t <- mm + nn
        # in total dim = mm + nn + 1
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [0.] * nn + [self.mu]
        
        # build G
        G_icol = list(range(mm)) * nn + list(range(mm, mm + nn)) \
            + [mm + nn] * nn 
        G_irow = [k  for k in range(nn) for _ in range(mm)] \
            + list(range(nn)) * 2 
        G_data = list(np.ravel(-self.rrate)) + [-1.] * nn \
            + [self.mu] * nn 
        G_icol += list(range(mm + nn + 1))   
        G_irow += list(range(nn, nn + mm + nn + 1))
        G_data += [-1.] * (mm + nn + 1)
        
        G_shape = (nn + mm + nn + 1, mm + nn + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn + mm + nn + 1)
        
        #build A
        A_icol = list(range(mm, mm + nn)) + list(range(mm)) + [mm + nn]
        A_irow = [0] * nn + [1] * (mm + 1)
        A_data = [1. / nn] * nn + [1.] * mm + [-1.]
        A_shape = (2, mm + nn + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # Omega
        self.sharpe = -res['pcost']
        # risk
        self.risk = 1. / res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][:mm] * self.risk)
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost'] * self.risk + self.mu
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
        
    def _sharpe_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # s <- [mm : mm + nn]
        # t <- mm + nn
        # in total dim = mm + nn + 1
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = [0.] * mm + [1. / nn] * nn + [0.]
        
        # build G
        G_icol = list(range(mm)) * nn + list(range(mm, mm + nn)) \
            + [mm + nn] * nn 
        G_irow = [k  for k in range(nn) for _ in range(mm)] \
            + list(range(nn)) * 2 
        G_data = list(np.ravel(-self.rrate)) + [-1.] * nn \
            + [self.mu] * nn 
        G_icol += list(range(mm + nn + 1))   
        G_irow += list(range(nn, nn + mm + nn + 1))
        G_data += [-1.] * (mm + nn + 1)
        
        G_shape = (nn + mm + nn + 1, mm + nn + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn + mm + nn + 1)
        
        #build A
        A_icol = list(range(mm)) + [mm + nn] + list(range(mm)) + [mm + nn]
        A_irow = [0] * (mm + 1) + [1] * (mm + 1)
        A_data = list(self.muk) +[-self.mu] + [1.] * mm + [-1.]
        A_shape = (2, mm + nn + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Omega
        self.sharpe = 1. / res['pcost']
        # risk
        self.risk =  res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _rr_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # s <- [mm : mm + nn]
        # in total dim = mm + nn
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [0.] * nn 
        
        # build G
        G_icol = list(range(mm)) * nn + list(range(mm, mm + nn))
        G_irow = [k  for k in range(nn) for _ in range(mm)] + list(range(nn)) 
        G_data = list(np.ravel(-self.rrate)) + [-1.] * nn 
        G_icol += list(range(mm + nn))   
        G_irow += list(range(nn, nn + mm + nn))
        G_data += [-1.] * (mm + nn)
        
        G_shape = (nn + mm + nn, mm + nn)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [-self.alpha[0]] * nn + [0.] * (mm + nn)
        
        #build A
        A_icol = list(range(mm, mm + nn)) + list(range(mm))
        A_irow = [0] * nn + [1] * mm
        A_data = [1. / nn] * nn + [1.] * mm 
        A_shape = (2, mm + nn)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
 
        # build b
        b_data = [self.risk, 1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # rate of return
        self.RR = -res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _risk_averse(self):
        # Order of variables
        # w <- [0 : mm]
        # s <- [mm : mm + nn]
        # in total dim = mm + nn
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [self.Lambda / nn] * nn
        
        # build G
        G_icol = list(range(mm)) * nn + list(range(mm, mm + nn)) 
        G_irow = [k  for k in range(nn) for _ in range(mm)] + list(range(nn)) 
        G_data = list(np.ravel(-self.rrate)) + [-1.] * nn 
        G_icol += list(range(mm + nn))   
        G_irow += list(range(nn, nn + mm + nn))
        G_data += [-1.] * (mm + nn)
        
        G_shape = (nn + mm + nn, mm + nn)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
  
        # build h
        h_data = [-self.alpha[0]] * nn + [0.] * (mm + nn)
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + nn)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # delta-risk
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
        