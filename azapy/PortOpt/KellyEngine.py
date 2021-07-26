# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 10:46:40 2021

@author: mircea
"""
import numpy as np
import pandas as pd
import cvxopt as cx
import warnings

class KellyEngine():
    """
    Computes the Kelly optimal portfolio.
    """
    def getWeights(self, rrate, rtype='Full', method='glpk'):
        """
        Computes the Kelly optimal weights.

        Parameters
        ----------
        rrate : pd.DataFrame
            Portfolio components historical rates of returns in the format 
           "date", "symbol1", "symbol2", etc. 
        rtype : string, optional
            Optimization type. it can be:\n
            'Full' - non-linear original Kelly problem. \n
            'Order2' - second order Taylor approximation of original Kelly 
            problem. It is a QP problem. \n
            The default is 'Full'.
        method : string, optional
            The QP solver class. It is relevant only if rtype='Order2'.
            It takes 2 values: 'glpk' or None for default cvxopt.solvers.qp 
            algorithm.
            The default is 'glpk'.

        Returns
        -------
        pd.Series
            Portfolio weights.
        """
        self.rrate = rrate
        self.method = method
        
        if rtype == 'Full':
            return self._calc_full()
        elif rtype == 'Order2':
            return self._calc_order2()
        else:
            raise ValueError("rtype must be either 'Full' or 'Order2")
            
    def _calc_full(self):
        nn = len(self.rrate.columns)
        
        # local function
        def F(x=None, z=None):
            if x is None: 
                return 0, cx.matrix(0.5/nn, (nn, 1))
            
            xx = np.array(x.T)[0]
        
            ff = self.rrate.apply(lambda row: np.dot(xx, row) + 1 , axis=1)
            val = -np.log(ff).mean()
            
            rf = self.rrate.apply(lambda col: col / ff)
            DF = cx.matrix(-rf.mean()).T
            
            if z is None: 
                return val, DF
            
            H = np.zeros((nn, nn))
            for i in range(nn):
                for j in range(i + 1):
                    H[i,j] = (rf.iloc[:,i] * rf.iloc[:,j]).mean() * z[0]
                    H[j,i] = H[i,j]
                    
            return val, DF, cx.matrix(H)
        
        icol = list(range(nn)) * 2
        irow = list(range(nn)) + [nn] * nn
        data = [-1.] * nn + [1.] * nn
        G = cx.spmatrix(data, irow, icol, (nn + 1, nn))
        h = cx.matrix([0.] * nn + [1.])
        dims = {'l': nn + 1, 'q': [], 's': []}
        
        res = cx.solvers.cp(F, G=G, h=h, dims=dims, 
                            options={'show_progress': False})
        
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        return pd.Series(res['x'], index=self.rrate.columns)
        
    def _calc_order2(self):
        mu = self.rrate.mean()
        nn = len(self.rrate.columns)
        rrd = self.rrate - mu

        # build P and q
        P = cx.matrix(0., (nn, nn))
        for i in range(nn):
            for j in range(i+1):
                P[i,j] = (rrd.iloc[:,i] * rrd.iloc[:,j]).mean()
                P[j,i] = P[i,j]
        
        q = cx.matrix(-mu)    

        #build G and h
        icol = list(range(nn)) * 2
        irow = list(range(nn)) + [nn] * nn
        data = [-1.] * nn + [1.] * nn
        G = cx.spmatrix(data, irow, icol, (nn + 1, nn))
        h = cx.matrix([0.] * nn + [1.])

        # solve
        res = cx.solvers.qp(P=P, q=q, G=G, h=h, solver=self.method, 
                            options={'show_progress': False})

        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            return pd.Series(np.nan, index=self.rrate.columns)      

        return pd.Series(res['x'], index=rrd.columns)
        