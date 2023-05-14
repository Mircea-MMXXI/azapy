from ._Port_Generator import _Port_Generator
from azapy.Generators.ModelPipeline import ModelPipeline
from azapy.Engines.KellyEngine import KellyEngine


class Port_Kelly(_Port_Generator):
    """
    Backtesting Kelly portfolio periodically rebalanced.
    
    **Attributes**
        * `pname` : `str` - portfolio name
        * `ww` : `pandasDataFrame` - portfolio weights at each rebalancing date
        * `port` : `pandas.Series` - portfolio historical time-series
        * `schedule` : `pandas.DataFrame` - rebalancing schedule
       
    The most important method is `set_model`. It must be called before any
    other method.
    """                  
    def set_model(self, rtype='ExpCone', hlength=1.25, method='ecos', 
                  verbose=False):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        rtype : `str`, optional
            Type of optimization. It could take the values:\n
                `'ExpCone'` - Exponential cone constraint programming solution
                for full Kelly problem. \n
                `'Full'` - Non-linear solver for full Kelly problem. \n
                `'Order2'` - Second order Tayler approximation of Kelly problem. \n
            The default is `'ExpCone'`.
        hlength : `float`, optional
            The length in year of the historical calibration period relative 
            to `'Dfix'`. A fractional number will be rounded to an integer number 
            of months. The default is `1.25` years. 
        method : `str`, optional
            The QP solver class. It is relevant only if `rtype='Order2'`.
            It takes 2 values: `'ecos'` or `'cvxopt'`.
            The default is `'ecos'`.
        verbose : Boolean, optional
            Sets verbose mode. The default is `False`.

        Returns
        -------
        `pandas.DataFrame` : The portfolio time-series in the format 'date', 
        'pcolname'.
        """
        mod = KellyEngine(colname=self.col_calib, freq=self.freq,
                          hlength=hlength, rtype=rtype, method=method)
        return super().set_model(ModelPipeline([mod]), verbose)