from ._Port_Generator import _Port_Generator
from azapy.Generators.ModelPipeline import ModelPipeline
from azapy.Analyzers.MVAnalyzer import MVAnalyzer


class Port_MV(_Port_Generator):
    """
    Backtesting MV (mean variance) portfolio periodically rebalanced.
    
    **Attributes**
        * `pname` : `str` - portfolio name
        * `ww` : `pandasDataFrame` - portfolio weights at each rebalancing date
        * `port` : `pandas.Series` - portfolio historical time-series
        * `schedule` : `pandas.DataFrame` - rebalancing schedule
      
    The most important method is `set_model`. It must be called before any
    other method.
    """                             
    def set_model(self, rtype='Sharpe', mu=None, mu0=0, aversion=None, 
                  ww0=None, hlength=3.25, method='ecos', verbose=False):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------.
        rtype : `str`, optional
            Optimization type. Possible values: \n
                `'Risk'` : optimal-risk portfolio for targeted expected rate of 
                return.\n
                `'Sharpe'` : Sharpe-optimal portfolio - maximization solution.\n
                `'Sharpe2'` : Sharpe-optimal portfolio - minimization solution.\n
                `'MinRisk'` : minimum risk portfolio.\n
                `'RiskAverse'` : optimal-risk portfolio for a fixed 
                risk-aversion factor.\n
                `'InvNrisk'` : optimal-risk portfolio with the same risk value 
                as a benchmark portfolio (e.g., same as equal weighted 
                portfolio).\n
                `'Diverse'` : optimal-diversified portfolio for targeted
                expected rate of return (maximum of inverse of 1-D).\n
                `'Diverse2'` : optimal-diversified portfolio for targeted
                expected rate of return (minimum of 1-D).\n
                `'MaxDiverse'` : maximum diversified portfolio.\n
                `'InvNdiverse'` : optimal-diversified portfolio with the same
                diversification factor as a benchmark portfolio 
                (e.g., same as equal weighted portfolio).\n
                `'InvNdrr'` : optima- diversified portfolio with the same 
                expected rate of return as a benchmark portfolio
                (e.g., same as equal weighted portfolio).\n
            The default is `'Sharpe'`.
        mu : `float`, optional
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype='Risk'`
            The default is `None`.
        mu0 : `float`, optional
            Risk-free rate accessible to the investor.
            Relevant only if `rtype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        aversion : `float`, optional
            The value of the risk-aversion factor.
            Must be positive. Relevant only if `rtype='RiskAvers'`.
            The default is `None`.
        ww0 : `list` (also `numpy.array` or `pandas.Series`), optional
            Targeted portfolio weights. 
            Relevant only if `rtype='InvNrisk'`.
            Its length must be equal to the number of
            symbols in `rrate` (`mktdata`). 
            All weights must be >= 0 with their sum > 0.
            If it is a list or a `numpy.array` then the weights are assumed to
            by in order of `rrate.columns`. If it is a `pandas.Series` then 
            the index should be compatible with the `rrate.columns` or `mktdata` 
            symbols (same symbols, not necessarily in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.
        hlength : `float`, optional
            The length in year of the historical calibration period relative
            to `'Dfix'`. A fractional number will be rounded to an integer 
            number of months. The default is `3.25` years.
        method : `str`, optional
            Quadratic programming numerical method. Could be `'ecos'` or
            `'cvxopt'`. The default is `'ecos'`.
        verbose : `Boolean`, optional
            If it set to `True` then it will print messages when the optimal
            portfolio degenerates to a single asset portfolio as a limited 
            case. 
            The default is `False`.

        Returns
        -------
        `pandas.DataFrame` : The portfolio time-series in the format 'date', 
        'pcolname'.
        """
        mod = MVAnalyzer(colname=self.col_calib, freq=self.freq,
                         hlength=hlength, rtype=rtype, mu=mu, d=1, mu0=mu0,
                         aversion=aversion, ww0=ww0, method=method)
        return super().set_model(ModelPipeline([mod]), verbose)