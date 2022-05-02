import numpy as np
import scipy.sparse as sps
import warnings

from .MADAnalyzer import MADAnalyzer
from ._solvers import _socp_solver

class LSSDAnalyzer(MADAnalyzer):
    """
    LSSD dispersion measure based portfolio optimization.
    
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
    def __init__(self, coef=[1.], 
                 mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, calendar=None,
                 rtype='Sharpe', method='ecos'):
        """
        Constructor

        Parameters
        ----------
        coef : list, optional
            List of coefficients (the list size defines the MAD 
            order). The default is [1.].
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
            SOCP numerical method. 
            Could be: 'ecos' or 'cvxopt'.
            The defualt is 'ecos'.

        Returns
        -------
        The object.
        """
        super().__init__(coef, mktdata, colname, freq, hlength, calendar,
                         rtype, method)
        
        
    def _set_method(self, method):
        self._set_socp_method(method)
        
        
    def _risk_calc_(self, prate):
        # mu = np.mean(prate)
        # prate = prate - mu
        # nn = len(prate)
        delta = []
        mu = 0.
        for _ in range(self.ll):
            xrate = mu - prate
            xrate[xrate <= 0.] = 0.
            dd = np.sqrt(np.mean(xrate**2))
            delta.append(dd)
            mu -= dd
            
        self.primary_risk_comp = np.array(delta)
        self.secondary_risk_comp = np.cumsum(self.primary_risk_comp) \
            -self.primary_risk_comp[0]
        self.risk = np.dot(self.primary_risk_comp, self.coef)
        
        return self.risk


    def _risk_min(self, d=1):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        nn = self.nn
        mm = self.mm
        ll = self.ll 
    
        # build c
        c_data = [0.] * mm 
        for l in range(ll):
            c_data += [self.coef[l]] + [0.] * nn
        
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += [k for k in range(l * nn, (l + 1) * nn) \
                       for _ in range(l)] \
                  + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * ((l + 1) * nn)
        G_icol += list(range(mm)) 
        G_irow += [ll * nn] * mm
        G_data += list(-self.muk * d)
        G_icol += list(range(mm))
        G_irow += list(range(nn * ll + 1, mm + nn * ll + 1))
        G_data += [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1) + 1, \
                                 mm + (l + 1) * (nn + 1)))
            G_irow += list(range(mm + ll * nn + 1 + l * nn, \
                                 mm + ll * nn + 1 + (l + 1) * nn))
            G_data += [-1.] * nn
        # cone
        for l in range(ll):
            G_icol += list(range((nn + 1) * l + mm , 
                                 (nn + 1) * l + mm + 1 + nn))
            G_irow += list(range(mm + 2 * nn * ll + 1 + l * (nn + 1), 
                                 mm + 2 * nn * ll + 1 + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (3 * nn * ll + 1 + mm + ll, mm + ll * (nn + 1))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        #build h
        h_data = [0.] * (nn * ll) + [-self.mu * d] \
               + [0.] * (mm + 2 * nn * ll + ll) 
        
        # build dims
        dims = {'l': (2 * ll * nn + 1 + mm), 'q': [nn + 1] * ll}
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # mLSSD
        self.risk = res['pcost']
        # LSSD
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primary_risk_comp) \
            - self.primary_risk_comp[0]
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        
        return self.ww
        
    
    def _sharpe_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        # u_l <- mm + l(nn+1), 
        # s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # and last t <- [mm + ll(nn + 1)]
        # in total dim = mm + ll(nn + 1) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [0.] * (ll * (nn + 1)) + [self.mu]
        
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += [k for k in range(l * nn, (l + 1) * nn) \
                       for _ in range(l)] \
                  + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * ((l + 1) * nn)
        G_icol += list(range(mm))
        G_irow += list(range(nn * ll, mm + nn * ll))
        G_data += [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1) + 1, \
                                 mm + (l + 1) * (nn + 1)))
            G_irow += list(range(mm + ll * nn + l * nn, \
                                 mm + ll * nn + (l + 1) * nn))
            G_data += [-1.] * nn
        # cone
        for l in range(ll):
            G_icol += list(range((nn + 1) * l + mm , 
                                 (nn + 1) * l + mm + 1 + nn))
            G_irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                                 mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (3 * nn * ll + mm + ll, mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (3 * nn * ll + mm + ll)
        
        # build dims
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll}
        
        # build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        A_icol += [mm + l * (nn + 1) for l in range(ll)]
        A_irow += [1] * ll
        A_data += list(self.coef)
        
        A_shape = (2, mm + ll * (nn + 1) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
 
        # build b
        b_data = [0., 1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # risk (=1/t)
        self.risk = 1. / res['x'][-1]
        # LSSD (=u)
        self.primary_risk_comp = np.array(
            [res['x'][mm + l * (nn + 1)] * self.risk for l in range(ll)])
        # Sharpe
        self.sharpe = -res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm] * self.risk)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
         # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primary_risk_comp) \
            - self.primary_risk_comp[0]
        
        return self.ww
    
    
    def _sharpe_inv_min(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        # u_l <- mm + l(nn+1), 
        # s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # and last t <- [mm + ll(nn + 1)]
        # in total dim = mm + ll(nn + 1) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = [0.] * mm 
        for l in range(ll):
            c_data += [self.coef[l]] + [0.] * nn
        c_data += [0.]
        
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += [k for k in range(l * nn, (l + 1) * nn) \
                       for _ in range(l)]\
                  + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * ((l + 1) * nn)
        G_icol += list(range(mm))
        G_irow += list(range(nn * ll, mm + nn * ll))
        G_data += [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1) + 1, \
                                 mm + (l + 1) * (nn + 1)))
            G_irow += list(range(mm + ll * nn + l * nn, \
                                 mm + ll * nn + (l + 1) * nn))
            G_data += [-1.] * nn
        # cone
        for l in range(ll):
            G_icol += list(range((nn + 1) * l + mm , 
                                 (nn + 1) * l + mm + 1 + nn))
            G_irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                                 mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (3 * nn * ll + mm + ll, mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (3 * nn * ll + mm + ll)
        
        # build dims
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll}
        
        # build A
        A_icol = (list(range(mm)) + [mm + ll * (nn + 1)]) * 2
        A_irow = [0] * (mm + 1) + [1] * (mm + 1)
        A_data = [1.] * mm + [-1.] + list(self.muk) + [-self.mu]
        A_shape = (2, mm + ll * (nn + 1) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
 
        # bjuild b
        b_data = [0., 1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # Sharpe
        self.sharpe = 1. / res['pcost']
        # risk 
        t = res['x'][-1]
        self.risk = res['pcost'] / t
        # LSSD (=u)
        self.primary_risk_comp = np.array(
            [res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primary_risk_comp) \
            - self.primary_risk_comp[0]
        
        return self.ww
        
    
    def _rr_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        nn = self.nn
        mm = self.mm
        ll = self.ll 
        
        # build c
        c_data = list(-self.muk) + [0.] * (ll * (nn + 1))
        
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += [k for k in range(l * nn, (l + 1) * nn) \
                       for _ in range(l)] \
                  + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * ((l + 1) * nn)
        G_icol += list(range(mm))
        G_irow += list(range(nn * ll, mm + nn * ll))
        G_data += [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1) + 1, \
                                 mm + (l + 1) * (nn + 1)))
            G_irow += list(range(mm + ll * nn + l * nn, \
                                 mm + ll * nn + (l + 1) * nn))
            G_data += [-1.] * nn
        # cone
        for l in range(ll):
            G_icol += list(range((nn + 1) * l + mm , 
                                 (nn + 1) * l + mm + 1 + nn))
            G_irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                                 mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (3 * nn * ll + mm + ll, mm + ll * (nn + 1))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (3 * nn * ll + mm + ll)
        
        # build dims
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll}
        
        # build A
        A_icol = list(range(mm)) + [mm + l * (nn + 1) for l in range(ll)]
        A_irow = [0] * mm + [1] * ll
        A_data = [1.] * mm + list(self.coef)
        A_shape = (2, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., self.risk]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # rate of return
        self.RR = -res['pcost']
        # LSSD
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 2)] \
                                           for l in range(ll)])
        # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primary_risk_comp) \
            - self.primary_risk_comp[0]
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm

        return self.ww
        
    
    def _risk_averse(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        nn = self.nn
        mm = self.mm
        ll = self.ll 
        
        # build c
        c_data = list(-self.muk)
        for l in range(ll):
            c_data += [self.Lambda * self.coef[l]] + [0.] * nn
        
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += [k for k in range(l * nn, (l + 1) * nn) \
                       for _ in range(l)]\
                  + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * ((l + 1) * nn)
        G_icol += list(range(mm))
        G_irow += list(range(nn * ll, mm + nn * ll))
        G_data += [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1) + 1, \
                                 mm + (l + 1) * (nn + 1)))
            G_irow += list(range(mm + ll * nn + l * nn, \
                                 mm + ll * nn + (l + 1) * nn))
            G_data += [-1.] * nn
        # cone
        for l in range(ll):
            G_icol += list(range((nn + 1) * l + mm , 
                                 (nn + 1) * l + mm + 1 + nn))
            G_irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                                 mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (3 * nn * ll + mm + ll, mm + ll * (nn + 1))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (3 * nn * ll + mm + ll)
        
        # build dims
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll}
        
        # build A 
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # mLSSD
        self.risk = (self.RR + res['pcost']) / self.Lambda
        # LSSD
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primary_risk_comp) \
            - self.primary_risk_comp[0]
        
        return self.ww
    