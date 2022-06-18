import numpy as np
import pandas as pd
import cvxopt as cx
import warnings

def gamblingKelly(pp=[0.6]):
    """
    Compute the Kelly allocation for multiple binary games.

    Parameters
    ----------
    `pp` : list, optional
        List of wining probabilities for each game. The default is `[0.6]`.

    Returns
    -------
    `pandas.Series`
        Bet sizes as percentage of the capital=1.
    """
    pp = np.array(pp)
    assert all((pp > 0) & (pp < 1)), "all pp must by in (0,1)"
    nn = len(pp)
    
    # local function 
    def _ffval(w, p):
        n=len(p)
        z = 0.
        df = np.zeros(n)
        dh = np.zeros((n,n))
        for x in range(2**n):
            ppro = 1.
            wsum = 1.
            for k in range(n):
                if x & (1 << k):
                    ppro *= p[k]
                    wsum += w[k]
                else:
                    ppro *= 1 - p[k]
                    wsum -= w[k]
                    
            z += ppro * np.log(wsum)
            for j in range(n):
                jsi = 1 if x & (1 << j) else -1
                term = jsi * ppro / wsum
                df[j] +=  term
                for i in range(j + 1):
                    isi = 1 if x & (1 << i) else -1
                    dh[i,j] -= isi * term / wsum
                    dh[j,i] = dh[i,j]
                    
        return z, df, dh
    
    # local function
    def _F(x=None, z=None):
        if x is None: 
            return 0, cx.matrix(0.5/nn, (nn, 1))
        
        xx = np.array(x.T)[0]
        # ensure that the log can be computed
        if np.sum(xx) >= 1.:
            xx = xx / np.sum(xx) * 0.999999
        f, df, dh = _ffval(xx, pp)
        val = -f
        DF = cx.matrix(-df).T
        if z is None: 
            return val, DF
        
        H = cx.matrix(-z[0] * dh)
        
        return val, DF, H

    # build G and h
    icol = list(range(nn)) * 2
    irow = list(range(nn)) + [nn] * nn
    data = [-1.] * nn + [1.] * nn
    G = cx.spmatrix(data, irow, icol, (nn + 1, nn))
    h = cx.matrix([0.] * nn + [1.])
    dims = {'l': nn + 1, 'q': [], 's': []}
    
    sol = cx.solvers.cp(F=_F, G=G, h=h, dims=dims, 
                        options={'show_progress': False})
    
    if 'optimal' not in sol['status']:
        warnings.warn(f"cannot find a good solution msg: {sol['status']}")
        
    return pd.Series(sol['x']).round(4) * 100
