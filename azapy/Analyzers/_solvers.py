"""
Wrappers for various solvers.
"""

import numpy as np
import scipy.sparse as sps
import scipy.optimize as sopt
import scipy.linalg as la
import warnings
import cvxopt as cx
import ecos

_tol_cholesky = 1.e-10

def _lp_scipy(c, G, h, A, b, method):
    # compute - suppress warning {'sparse': True}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = sopt.linprog(c=c, A_ub=G, b_ub=h, 
                           A_eq=A, b_eq=b, method=method, 
                           options={'sparse': True})
        
    # gather the results
    rout = {}
    rout['status'] = res.status
    rout['infostring'] = res.message
    rout['pcost'] = np.float64(res.fun)
    rout['x'] = res.x
    
    return rout

def _lp_cvxopt(c, G, h, A, b, method):
    c = cx.matrix(c)
    G = cx.spmatrix(G.data, G.row, G.col, size=G.shape)
    h = cx.matrix(h)
    if A is not None:
        A = cx.spmatrix(A.data, A.row, A.col, size=A.shape)
        b = cx.matrix(b)
    
    if method == 'cvxopt': 
        solver = None
    else: 
        solver = method
    
    res = cx.solvers.lp(c, G, h, A, b, solver=solver, 
                            options={'show_progress': False,
                                     'glpk':{'msg_lev':'GLP_MSG_OFF'}})
    
    # gather the results
    rout = {}
    rout['status'] = 0 if 'optimal' in res['status'] else 2
    rout['infostring'] = res['status']
    rout['pcost'] = np.float64(res['primal objective'])
    rout['x'] = res['x']
    
    return rout

def _lp_ecos(c, G, h, A, b, method):
    c = np.array(c)
    G = G.tocsc()
    h = np.array(h)
    if A is not None:
        A = A.tocsc()
        b = np.array(b)
    dims = {'l': G.shape[0], 'q': [], 'e': 0}
    
    kwargs = {'verbose': False}
    res = ecos.solve(c, G, h, dims, A, b, **kwargs)
    
    # gather the results
    rout = {}
    rout['status'] = res['info']['exitFlag']
    # accept 10 Close to optimal as optimal
    if rout['status'] == 10: rout['status'] = 0
    rout['infostring'] = res['info']['infostring']
    rout['pcost'] = np.float64(res['info']['pcost'])
    rout['x'] = res['x']
    
    return rout
    
def _lp_solver(method, c_data, G, h_data, A=None, b_data=None):
    if method == 'ecos':
        return _lp_ecos(c_data, G, h_data, A, b_data, method)
    elif method in ['glpk', 'conelp', 'cvxopt']:
        return _lp_cvxopt(c_data, G, h_data, A, b_data, method)
    else:
        return _lp_scipy(c_data, G, h_data, A, b_data, method)

def _socp_cvxopt(c, G, h, dims, A, b):
    c = cx.matrix(c)
    G = cx.spmatrix(G.data, G.row, G.col, size=G.shape)
    h = cx.matrix(h)
    dims['s'] = []
    if A is not None:
        A = cx.spmatrix(A.data, A.row, A.col, size=A.shape)
        b = cx.matrix(b)
        
    res = cx.solvers.conelp(c, G, h, dims, A, b, 
                            options={'show_progress': False})
    
    # gather the results
    rout = {}
    rout['status'] = 0 if 'optimal' in res['status'] else 2
    rout['infostring'] = res['status']
    rout['pcost'] = np.float64(res['primal objective'])
    rout['x'] = res['x']
    
    return rout

def _socp_ecos(c, G, h, dims, A, b):
    c = np.array(c)
    G = G.tocsc()
    h = np.array(h)
    dims['e'] = 0
    if A is not None:
        A = A.tocsc()
        b = np.array(b)
    
    kwargs = {'verbose': False, 'max_iters': 1000}
    res = ecos.solve(c, G, h, dims, A, b, **kwargs)
    
    # gather the results
    rout = {}
    rout['status'] = res['info']['exitFlag']
    # accept 10 Close to optimal as optimal
    if rout['status'] == 10: rout['status'] = 0
    rout['infostring'] = res['info']['infostring']
    rout['pcost'] = np.float(res['info']['pcost'])
    rout['x'] = res['x']
    
    return rout

def _socp_solver(method, c_data, G, h_data,dims, A=None, b_data=None):
    if method == 'ecos':
        return _socp_ecos(c_data, G, h_data, dims, A, b_data)
    else:
        return _socp_cvxopt(c_data, G, h_data, dims, A, b_data)

def _qp_cvxopt(P, q, G, h, A, b):
    P = cx.matrix(P)
    q = cx.matrix(q)
    G = cx.spmatrix(G.data, G.row, G.col, size=G.shape)
    h = cx.matrix(h)
    if A is not None:
        A = cx.spmatrix(A.data, A.row, A.col, size=A.shape)
        b = cx.matrix(b)
    # else:
    #     AA = None
    #     bb = None
    
    res = cx.solvers.qp(P, q, G, h, A, b, 
                        options={'show_progress': False})
    
    # gather the results
    rout = {}
    rout['status'] = 0 if 'optimal' in res['status'] else 2
    rout['infostring'] = res['status']
    rout['pcost'] = np.float64(res['primal objective'])
    rout['x'] = res['x']
    
    return rout

def _qp_ecos(P, q, G, h, A, b):
    cc = np.concatenate((q, [1.]))
    
    gr, gc = G.shape
    pr, pc = P.shape
    
    gg = sps.block_diag((G, [-1.]))
    
    if any(np.diag(P) < _tol_cholesky):
        pp = sps.block_diag((-la.sqrtm(P), [1.]))
    else:
        pp = sps.block_diag((-la.cholesky(P, overwrite_a=True), [1.]))
        
    pp = sps.block_diag((-la.sqrtm(P), [1.]))
    GG = sps.vstack([gg, pp], format='csc')
  
    dims = {'l': gr, 'q': [pr + 2], 'e': 0}
    
    hh = np.concatenate((h, [1.] + [0.] * (pr + 1)))

    if A is not None:
        A = sps.hstack((A, np.zeros((A.shape[0], 1))), format='csc')
        b = np.array(b)

    kwargs = {'verbose': False}
    res = ecos.solve(cc, GG, hh, dims, A, b, **kwargs)
    
    # gather the results
    rout = {}
    rout['status'] = res['info']['exitFlag']
    # accept 10 Close to optimal as optimal
    if rout['status'] == 10: rout['status'] = 0
    rout['infostring'] = res['info']['infostring']
    rout['pcost'] = np.float64(res['info']['pcost'] + 0.5)
    rout['x'] = res['x'][:-1]
    
    return rout
    
def _qp_solver(method, P, q_data, G, h_data, A=None, b=None):
    if method == 'ecos':
        return _qp_ecos(P, q_data, G, h_data, A, b)
    else:
        return _qp_cvxopt(P, q_data, G, h_data, A, b)
    
    
def _exp_cone_ecos(c, G, h, dims, A, b):
    c = np.array(c)
    G = G.tocsc()
    h = np.array(h)

    if A is not None:
        A = A.tocsc()
        b = np.array(b)
    
    kwargs = {'verbose': False, 'max_iters': 1000}
    res = ecos.solve(c, G, h, dims, A, b, **kwargs)
    
    # gather the results
    rout = {}
    rout['status'] = res['info']['exitFlag']
    # accept 10 Close to optimal as optimal
    if rout['status'] == 10: rout['status'] = 0
    rout['infostring'] = res['info']
    rout['pcost'] = np.float64(res['info']['pcost'])
    rout['x'] = res['x']
 
    return rout


def _exp_cone_solver(method, c, G, h, dims, A=None, b=None):
    # only ecos so far - don't check for method
    return _exp_cone_ecos(c, G, h, dims, A, b)
    