import numpy as np
import scipy.sparse as sps
import warnings

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _lp_solver

class MADAnalyzer(_RiskAnalyzer):
    """
    m-level MAD based optimal portfolio strategies.

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
    def __init__(self, coef=[1.], mktdata=None, colname='adjusted', freq='Q',
                 hlength=3.25, calendar=None, rtype='Sharpe', method='ecos'):
        """
        Constructor

        Parameters
        ----------
        `coef` : list, optional
            Positive non-increasing list of mixture coefficients. 
            The default is `[1.]`.
        `mktdata` : `pandas.DataFrame`, optional
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        `colname` : str, optional
            Name of the price column from mktdata used in the weights
            calibration. The default is `'adjusted'`.
        `freq` : str, optional
            Rate of return horizon in number of business day. it could be
            'Q' for quarter or 'M' for month. The default is `'Q'`.
        `hlength` : float, optional
            History length in number of years used for calibration. A
            fractional number will be rounded to an integer number of months.
            The default is `3.25` years.
        `calendar` : `numpy.busdaycalendar`, optional
            Business days calendar. If is it `None` then the calendar will
            be set to NYSE business calendar.
            The default is `None`.
        `rtype` : str, optional
            Optimization type. Possible values \n
                'Risk' : minimization of dispersion (risk) measure for 
                targeted rate of return. \n
                'Sharpe' : maximization of generalized Sharpe ratio.\n
                'Sharpe2' : minimization of the inverse generalized Sharpe 
                ratio.\n
                'MinRisk' : minimum dispersion (risk) portfolio.\n
                'InvNrisk' : optimal portfolio with the same dispersion (risk)
                value as a benchmark portfolio 
                (e.g. equal weighted portfolio).\n
                'RiskAverse' : optimal portfolio for a fixed value of 
                risk-aversion factor.
            The default is `'Sharpe'`.
        `method` : str, optional
            method : string, optional
            Linear programming numerical method.
            Could be:'ecos', 'highs-ds', 'highs-ipm', 'highs',
            'interior-point', 'glpk' and 'cvxopt'.
            The defualt is `'ecos'`.

        Returns
        -------
        The object.

        """
        super().__init__(mktdata, colname, freq, hlength, calendar, rtype)

        self._set_method(method)

        self.coef = np.array(coef)
        if any(self.coef <= 0.):
            raise ValueError("All coef must be positive")
        if any(i < j for i, j in zip(self.coef, self.coef[1:])):
            raise ValueError("coef list should not be increasing")
            
        self.ll = self.coef.size
        
        if self.ll <= 0:
            raise ValueError("coef must contain at least one element")
            
        self.coef = self.coef / np.sum(self.coef)
        
        
    def _set_method(self, method):
        self._set_lp_method(method)


    def getRisk(self, ww, rrate=None):
        """
        Returns the value of MAD for a give portfolio.

        Parameters
        ----------
        `ww` : list (`numpy.array` or `pandas.Series`)
            Portfolio weights.
        `rrate` : `pandas.series`, optional
            The portfolio components historical rates of returns.
            If it is not `None`, it will overwrite the rrate computed in the
            constructor from mktdata. The default is `None`.

        Returns
        -------
        float :
        The value of mMAD
        """
        if rrate is not None:
            self.set_rrate(rrate)

        w = np.array(ww)
        if any(w < 0.):
            raise ValueError("All ww must be non negative")
        w = w / w.sum()

        prate = np.dot(self.rrate, w)

        self._risk_calc_(prate)
        self.RR = np.dot(w, self.muk)
        self.ww = w

        return self.risk


    def _risk_calc_(self, prate):
        # mu = np.mean(prate)
        # prate = prate - mu
        # nn = len(prate)
 
        delta = []
        mu = 0.
        for _ in range(self.ll):
            xrate = mu - prate
            xrate[xrate <= 0.] = 0.
            dd = np.mean(xrate)
            delta.append(dd)
            mu -= dd

        self.primary_risk_comp = np.array(delta)
        self.secondary_risk_comp = \
            np.cumsum(np.insert(self.primary_risk_comp, 0, 0))[:-1]
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
        G_irow += [ll * nn] * mm
        G_data += list(-self.muk * d)
        
        G_icol += list(range(mm + ll * (nn + 1)))
        G_irow += list(range(ll * nn + 1, ll * nn + 1 + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1))

        G_shape = (nn * ll + 1 + mm + ll * (nn + 1), mm + (nn + 1) * ll)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)

        # build h
        h_data = [0.] * (nn * ll) + [-self.mu * d] \
               + [0.] * (mm + ll * (nn + 1))

        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        for l in range(ll):
            A_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            A_irow += [l + 1] * (nn + 1)
            A_data += [-1.] + [1./nn] * nn

        A_shape = (ll + 1, mm + (nn + 1) * ll)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1.] + [0.] * ll

        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)

        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)

        # mMAD
        self.risk = res['pcost']
        # delta-risk values
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # delta-risk strikes
        self.secondary_risk_comp = \
            np.cumsum(np.insert(self.primary_risk_comp, 0, 0))[:-1]
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
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) * (l + 1)
            G_data += [-1.] * ((l + 1) * nn)
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1) + 1))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)

        G_shape = (nn * ll + mm + ll * (nn + 1) + 1, mm + (nn + 1) * ll + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)

        # build h
        h_data = [0.] * (nn * ll + mm + ll * (nn + 1) + 1)

        # build A
        A_icol = [mm + l * (nn + 1) for l in range(ll)]
        A_irow = [0] * ll
        A_data = list(self.coef)
        for l in range(ll):
            A_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            A_irow += [l + 1] * (nn + 1)
            A_data += [-1.] + [1. / nn] * nn
        
        A_icol += list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow += [ll + 1] * (mm + 1)
        A_data += [1.] * mm + [-1.]

        A_shape = (ll + 2, mm + (nn + 1) * ll + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1.] + [0.] * (ll + 1)

        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)

        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)

        t = res['x'][-1]
        # mMAD-Sharpe
        self.sharpe = -res['pcost']
        # mMAD
        self.risk = 1 / t
        # delta-risk values
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] / t \
                                           for l in range(ll)])
        # delta-risk strikes
        self.secondary_risk_comp = \
            np.cumsum(np.insert(self.primary_risk_comp, 0, 0))[:-1]
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost'] / t + self.mu

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
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += \
                [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)] \
                + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * ((l + 1) * nn)
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1) + 1))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)

        G_shape = (nn * ll + mm + ll * (nn + 1) + 1, mm + (nn + 1) * ll + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)

        # build h
        h_data = [0.] * (nn * ll + mm + ll * (nn + 1) + 1)

        # build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow = [0] * (mm + 1)
        A_data = list(self.muk) + [-self.mu]
        for l in range(ll):
            A_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            A_irow += [l + 1] * (nn + 1)
            A_data += [-1.] + [1. / nn] * nn
        
        A_icol += list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow += [ll + 1] * (mm + 1)
        A_data += [1.] * mm + [-1.]

        A_shape = (ll + 2, mm + (nn + 1) * ll + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1.] + [0.] * (ll + 1)

        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)

        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)

        t = res['x'][-1]
        # mMAD-Sharpe
        self.sharpe = 1. / res['pcost']
        # mMAD
        self.risk = res['pcost'] / t
        # delta-risk values
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] / t \
                                           for l in range(ll)])
        # delta-risk strikes
        self.secondary_risk_comp = \
            np.cumsum(np.insert(self.primary_risk_comp, 0, 0))[:-1]
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu

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
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += \
                  [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)] \
                  + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * ((l + 1) * nn)
            
        G_icol += list(range(mm + ll * (nn + 1)))
        G_irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1))

        G_shape = (nn * ll + mm + ll * (nn + 1), mm + (nn + 1) * ll)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)

        # build h
        h_data = [0.] * (nn * ll + mm + ll * (nn + 1))

        # build A
        A_icol = [mm + l * (nn + 1) for l in range(ll)]
        A_irow = [0] * ll
        A_data = list(self.coef)
        
        A_icol += list(range(mm))
        A_irow += [1] * mm
        A_data += [1.] * mm
        for l in range(ll):
            A_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            A_irow += [l + 2] * (nn + 1)
            A_data += [-1.] + [1. / nn] * nn

        A_shape = (ll + 2, mm + (nn + 1) * ll)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [self.risk, 1.]  + [0.] * ll

        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)

        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)

        # delta-risk values
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # delta-risk strikes
        self.secondary_risk_comp = \
            np.cumsum(np.insert(self.primary_risk_comp, 0, 0))[:-1]
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost']

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
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += [k for k in range(l * nn, (l + 1) * nn) \
                       for _ in range(l)] + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * ((l + 1) * nn)
            
        G_icol += list(range(mm + ll * (nn + 1)))
        G_irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1))

        G_shape = (nn * ll + mm + ll * (nn + 1), mm + (nn + 1) * ll)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)

        # build h
        h_data = [0.] * (nn * ll + mm + ll * (nn + 1))

        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        for l in range(ll):
            A_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            A_irow += [l + 1] * (nn + 1)
            A_data += [-1.] + [1. / nn] * nn
        
        A_shape = (ll + 1, mm + (nn + 1) * ll)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1.] + [0.] * ll

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
        # mMAD
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # delta-risk values
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # delta-risk strikes
        self.secondary_risk_comp = \
            np.cumsum(np.insert(self.primary_risk_comp, 0, 0))[:-1]

        return self.ww
