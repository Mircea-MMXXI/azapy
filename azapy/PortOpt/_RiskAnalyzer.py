# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:29:33 2021

@author: mircea 
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

class _RiskAnalyzer:
    """
    Base class for all XXXAnalyzer classes.
    All derive classes need to implement:
        _risk_calc(self, prate, alpha)
        _risk_min(self, d=1)
        _sharpe_max(self)
        _sharpe_inv_min(self)
        _rr_max(self)
        _risk_averse(self)
    """
    
    def __init__(self, rrate, rtype):
        """
        Constructor

        Parameters
        ----------
        rrate : pandas.DataFrame, optional
            Portfolio components historical rates of returns in the format 
           "date", "symbol1", "symbol2", etc. The default is None.
        rtype : string, optional
            Optimization type. Possible values \n
                "Risk" : minimization of dispersion (risk) measure.\n
                "Sharpe" : maximization of generalized Sharpe ratio.\n
                "Sharpe2" : alternative computation of generalized Sharpe 
                ratio.\n
                "MinRisk" : optimal portfolio with minimum dispersion (risk) 
                value.\n
                "InvNRisk" : optimal portfolio with the same dispersion (risk)
                value as equally weighted portfolio. \n
                "RiskAverse" : optimal portfolio for a fixed risk aversion 
                coefficient.
            The default is "Sharpe".
        
        Returns
        -------
        The object.
        """
        self.ww = None
        self.status = None
        
        self.rrate = None
        self.nn = None
        self.mm = None
        self.muk = None
        self.mu = None
        self.alpha = [1.]
        self.coef = [1.]
        
        self.risk = None
        self.primery_risk_comp = None
        self.secondary_risk_comp = None
        self.sharpe = None
        self.RR = None
        
        self.rtype = None
        
        self.set_rtype(rtype)
        
        if rrate is not None:
            self.set_rrate(rrate)
            
        self.rng = None
        self.set_random_seed()
        
    def getWeights(self, mu, rrate=None, rtype=None, d=1):
        """
        Computes the optimal portfolio weights.

        Parameters
        ----------
        mu : float
            Rate of reference. Its meaning depends on the optimization 
            criterion. For rtype set to\n
                "Risk" : mu is the targeted portfolio rate of returns.\n
                "Sharpe" and "Sharpe2" : mu is the risk-free rate.\n
                "MinRisk" and "InvNRisk": mu is ignored. \n
                "RiskAverse" : mu is the Lambda aversion coefficient.
        rrate : pandas.DataFrame, optional
            MkT Data. If is not None it will overwrite the rrate set by the 
            constructor. The default is None.
        rtype : string, optional
            Optimization type. If is not None it will overwrite the value 
            set by the constructor. The default is None.
        d : int, optional
            Frontier type. Has effect only if rtype="Risk". A value of 1 will 
            trigger the evaluation of optimal portfolio along the efficient 
            frontier. Otherwise it will find the portfolio with the lowest 
            rate of return along the inefficient portfolio frontier. 
            The default is 1.
        
        Returns
        -------
        pandas.Series
            Portfolio weights. 
        """
        if rrate is not None:
            self.set_rrate(rrate)
            
        if rtype is not None:
            self.set_rtype(rtype)
        
        if self.rtype == "Risk":
            if mu > self.muk.max(): 
                self.mu = self.muk.max()
            elif mu < self.muk.min():
                self.mu = self.muk.min()
            else:
                self.mu = mu
 
            self._risk_min(d=d)
        elif self.rtype == "Sharpe":
            if mu > self.muk.max(): 
                self.mu = self.muk.max() * 0.999
            else:
                self.mu = mu
                
            self._sharpe_max()
        elif self.rtype == "MinRisk":
            self.mu = self.muk.min()
            
            self._risk_min(d=1)
        elif self.rtype == "Sharpe2":
            if mu > self.muk.max(): 
                self.mu = self.muk.max() * 0.999
            else:
                self.mu = mu
                
            self._sharpe_inv_min()
        elif self.rtype == "InvNrisk":
            ww = np.array([1.] * len(self.rrate.columns))
            ww = ww / np.sum(ww)
            self.getRisk(ww)
            self._rr_max()
        elif self.rtype == "RiskAverse":
            self.Lambda = mu
            self._risk_averse()
        else:
            print(f"should not be here!! Unknown rtype {rtype}")
            return np.nan
                     
        self.ww = pd.Series(self.ww, index=self.rrate.columns)
        return self.ww
    
    def getRisk(self, ww, rrate=None):
        """
        Returns the value of the dispersion (risk) measure for a give portfolio.

        Parameters
        ----------
        ww : list (np.array or pandas.Series)
            Portfolio weights. Must have a size equal to the number of symbols 
            in rrate (MkT data). All weights must by >= 0 with at least one 
            > 0.
        rrate : pandas.series, optional
            MkT Data. If is not None it will overwrite the rrate set by the 
            constructor. The default is None.

        Returns
        -------
        float
        The dispersion (risk) measure value.

        """
        if rrate is not None: 
            self.set_rrate(rrate)
            
        w = np.array(ww)
        assert all(w >= 0.), "All weights must be non-negative"
        w = w / w.sum()
        
        prate = np.dot(self.rrate, w)
        
        self.primery_risk_comp = []
        self.secondary_risk_comp = []
        for alpha in self.alpha:
            self.status, secondary, primery = self._risk_calc(prate, alpha)
            if self.status != 0:
                print(
                    f"status: {self.status}, wrong risk calc for alpha {alpha}"
                    )
                return np.nan
            
            self.primery_risk_comp.append(primery)
            self.secondary_risk_comp.append(secondary)
        self.risk = np.dot(self.primery_risk_comp, self.coef)
        self.primery_risk_comp = np.array(self.primery_risk_comp)
        self.secondary_risk_comp = np.array(self.secondary_risk_comp)
        self.RR = np.dot(w, self.muk)
        self.ww = w
        
        return self.risk
    
    def set_rrate(self, rrate):
        """
        Sets portfolio components historical rates of returns in the format 
        "date", "symbol1", "symbol2", etc. 

        Parameters
        ----------
        rrate : pandas.DataFrame
            MkT Data. It will overwrite the value set by the constructor.
        Returns
        -------
        None.
        """
        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean()
        self.rrate = rrate - self.muk
        
    def set_rtype(self, rtype):
        """
        Sets the optimization type.

        Parameters
        ----------
        rtype : string
            Optimization type. It will overwrite the value set by the 
            constructor.
        Returns
        -------
        None.
        """
        rtypes = ["Sharpe", "Risk", "MinRisk", "Sharpe2", "InvNrisk", 
                  "RiskAverse"]
        assert rtype in rtypes, f"type must be one of {rtypes}"
        self.rtype = rtype
        
    def viewFrontiers(self, efficient=20, inefficient=20, musharpe=0., 
                      component=True, randomport=20, inverseN=True,
                      fig_type='RR_risk',
                      options=None, save=None, data=None):
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
        component : boolean, optional
            If True the portfolios containing a single component are evaluated 
            and added to the plot for references. The default is True.
        randomport : int, optional
            The number of portfolios with random weights (inefficient) to be 
            evaluate and added to the plot for reference. The default is 20.
        inverseN : boolean, optional
            If True the equally weighted portfolio and the optimal portfolio 
            with the same dispersion (risk) value are evaluated and added to 
            the plot. The default is True.
        fig_type : string, optional
            Graphical representation format.
            If it is set to "RR_risk" the data is plotted in the rate of return 
            vs dispersion representation, otherwise the Sharpe vs rate of 
            return will be used. The default is 'RR_risk'.
        options : dictionary, optional
            Additional graphical setups (keys): "title", "xlabel", "ylabel", 
            "tangent".\n
            "title", "xlabel" and "ylabel" are strings overwriting the default 
            values. \n
            "tangent" is a boolean. If set to True it will print
            the Sharpe tangent. The default is True.
        save : string, optional
            File name to save the plot. The default is None.
        data : dictionary, optional
            Numerical data to construct the plot. If it is not None it 
            will take precedence and no other numerical evaluation will be 
            performed. It is meant to produce different plot representations
            without recomputation. The default is None.

        Returns
        -------
        dictionary
            Numerical data used to make the plots. It can be passed back to
            reconstruct the plots without reevaluations.
        """
        if data is not None:
            data['save'] = save
            if fig_type == 'RR_risk':
                self._plot_f1(data)
            else:
                self._plot_f2(data)
            return data
            
        res = defaultdict(lambda: None)
        res['options'] = options
        res['save'] = save
        # min risk
        res['risk_min'] = defaultdict(lambda x=0: None)
        rr = self.muk.min()
        ww_min = self.getWeights(mu=rr, rtype='Risk')
        risk_min = self.risk
        rr_min = self.RR
        res['risk_min']['risk'] = risk_min
        res['risk_min']['rr'] = rr_min
        res['risk_min']['ww'] = ww_min
        
        # efficient frontier
        res['efficient'] = defaultdict(lambda: None)
        res['efficient']['efficient'] = efficient
        if efficient > 0:
            rr_eff = np.linspace(rr_min, self.muk.max(), efficient)
            risk_eff = []
            ww_eff = []
            RR_eff = []
            for rr in rr_eff:
                ww_eff.append(self.getWeights(mu=rr, rtype='Risk'))
                risk_eff.append(self.risk)
                RR_eff.append(self.RR)
            res['efficient']['risk'] = np.array(risk_eff)
            res['efficient']['rr'] = np.array(RR_eff)
            res['efficient']['ww'] = ww_eff
                
        # inefficient frontier
        res['inefficient'] = defaultdict(lambda: None)
        res['inefficient']['inefficient'] = inefficient
        if inefficient > 0 :
            rr_ineff = np.linspace(self.muk.min(), rr_min, inefficient)
            risk_ineff = []
            ww_ineff = []
            RR_ineff = []
            for rr in rr_ineff:
                ww_ineff.append(self.getWeights(mu=rr, rtype='Risk', d=-1))
                risk_ineff.append(self.risk) 
                RR_ineff.append(self.RR)
            res['inefficient']['risk'] = np.array(risk_ineff)
            res['inefficient']['rr'] = np.array(RR_ineff)
            res['inefficient']['ww'] = ww_ineff
                
        # sharpe point
        res['sharpe'] = defaultdict(lambda: None)
        res['sharpe']['mu'] = musharpe
        if not np.isnan(musharpe):
            if musharpe > self.muk.max(): musharpe = self.muk.max() * 0.9999
            ww_sharpe = self.getWeights(mu=musharpe, rtype='Sharpe')
            rr_sharpe = self.RR
            risk_sharpe = self.risk
            sharpe = self.sharpe
            res['sharpe']['mu'] = musharpe
            res['sharpe']['risk'] = risk_sharpe
            res['sharpe']['rr'] = rr_sharpe
            res['sharpe']['ww'] = ww_sharpe
            res['sharpe']['sharpe'] = sharpe
            
        # portfolio components
        res['component'] = defaultdict(lambda: None)
        res['component']['component'] = component
        if component:
            risk_comp = []
            rr_comp = []
            label_comp = []
            ww_comp = []
            for k in range(self.mm):
                ww = [0.] * self.mm
                ww[k] = 1.
                risk_comp.append(self.getRisk(ww))
                rr_comp.append(self.muk[k])
                label_comp.append(self.muk.index[k])
                ww_comp.append(ww)
            res['component']['risk'] = np.array(risk_comp)
            res['component']['rr'] = np.array(rr_comp)
            res['component']['ww'] = ww_comp
            res['component']['label'] = label_comp
            
        # random portfolios
        res['randomport'] = defaultdict(lambda: None)
        res['randomport']['randomport'] = randomport
        if randomport > 0:
            risk_rp = []
            rr_rp = []
            ww_rp = []
            for _ in range(randomport):
                ww = self._ww_gen()
                risk_rp.append(self.getRisk(ww))
                rr_rp.append(self.RR)
                ww_rp.append(ww)
            res['randomport']['risk'] = np.array(risk_rp)
            res['randomport']['rr'] = np.array(rr_rp)
            res['randomport']['ww'] = ww_rp
            
        # inverse N
        res['inverseN'] = defaultdict(lambda: None)
        res['inverseN']['inverseN'] = inverseN
        if inverseN:
            risk_n = []
            rr_n = []
            ww_n = []
            # 1/N portfolio
            wwn = np.array([1.] * self.mm)
            wwn = wwn / np.sum(wwn)
            risk_n.append(self.getRisk(wwn))
            rr_n.append(self.RR)
            ww_n.append(wwn)
            # optimal 1/N risk portfolio
            ww_n.append(self.getWeights(0., rtype="InvNrisk"))
            risk_n.append(self.risk)
            rr_n.append(self.RR)
            res['inverseN']['risk'] = np.array(risk_n)
            res['inverseN']['rr'] = rr_n
            res['inverseN']['ww'] = ww_n
            res['inverseN']['label'] = ['1/N', 'InvNRisk']
            
        # Plot
        if fig_type == 'RR_risk':
            self._plot_f1(res)
        else:
            self._plot_f2(res)
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
            xl1, xl2 = plt.xlim()
            yl1, yl2 = plt.ylim()
            lff = 1.02
            xl2 *= lff
            yl2 *= lff
            ax.set(xlim=(xl1, xl2), ylim=(yl1, yl2))
            
        if not np.isnan(res['sharpe']['mu']):
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
   
        if res['randomport']['randomport'] > 0:
            ax.scatter(x=res['randomport']['risk'], y=res['randomport']['rr'],
                       marker='o', s=9, edgecolors='b', alpha=0.5 )
            
        if res['inverseN']['inverseN']:
            ax.scatter(x=res['inverseN']['risk'], y=res['inverseN']['rr'],
                      marker='x', c='g', s=36)
            lf = 1.01
            for k in range(len(res['inverseN']['label'])):
                ax.annotate(res['inverseN']['label'][k],
                            (res['inverseN']['risk'][k] * lf, 
                             res['inverseN']['rr'][k] * lf))
            
        if res['save'] is not None:
            plt.savefig(res['save'])
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
        
        musharpe = np.where(np.isnan(res['sharpe']['mu']), 0., 
                            res['sharpe']['mu'])
        
        shp = (res['risk_min']['rr'] - musharpe) / res['risk_min']['risk']
        ax.scatter(x=res['risk_min']['rr'], y=shp, marker='D', c='g', s=36)
        
        if res['efficient']['efficient'] > 0:
            shp = (res['efficient']['rr'] - musharpe) \
                / res['efficient']['risk']
            ax.plot(res['efficient']['rr'], shp,  
                    color='lightblue', linewidth=2)
            
        if res['inefficient']['inefficient'] > 0:
            shp = (res['inefficient']['rr'] - musharpe) \
                / res['inefficient']['risk']
            ax.plot(res['inefficient']['rr'], shp,
                    color='red', linewidth=2)
            
        if res['component']['component']:
            lf = 1.01
            for k in range(len(res['component']['risk'])):
                shp = (res['component']['rr'][k] - musharpe) \
                    / res['component']['risk'][k]
                ax.scatter(x=res['component']['rr'][k], y=shp, 
                           marker='s', c='b', s=16)
                ax.annotate(res['component']['label'][k],
                            (res['component']['rr'][k] * lf, shp * lf))
            xl1, xl2 = plt.xlim()
            yl1, yl2 = plt.ylim()
            lff = 1.02
            xl2 *= lff
            yl2 *= lff
            ax.set(xlim=(xl1, xl2), ylim=(yl1, yl2))
            
        if not np.isnan(res['sharpe']['mu']):
            ax.scatter(x=res['sharpe']['rr'], y=res['sharpe']['sharpe'], 
                       marker='D', c='g', s=26)
            
        if res['randomport']['randomport'] > 0:
            shp = (res['randomport']['rr'] - musharpe) \
                / res['randomport']['risk']
            ax.scatter(x=res['randomport']['rr'], y=shp,
                       marker='o', s=9, edgecolors='b', alpha=0.5 )
            
        if res['inverseN']['inverseN']:
            shp = (res['inverseN']['rr'] - musharpe) \
                / res['inverseN']['risk']
            ax.scatter(x=res['inverseN']['rr'], y=shp,
                      marker='x', c='g', s=36)
            lf = 1.01
            for k in range(len(res['inverseN']['label'])):
                ax.annotate(res['inverseN']['label'][k],
                            (res['inverseN']['rr'][k] * lf, shp[k] * lf))
            
        if res['save'] is not None:
            plt.savefig(res['save'])
        plt.show()
        
    def set_random_seed(self, seed = 42):
        """
        Sets the seed for Dirichlet random generator used in viewFrontiers.

        Parameters
        ----------
        seed : int, optional
            The random generator seed in case you want to set it to a weird 
            value other than 42 :). The default is 42.

        Returns
        -------
        None.
        """
        self.rng = np.random.RandomState(seed)
        
    def _ww_gen(self):
        return self.rng.dirichlet([0.5] * self.mm)
        
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
        