from ._Port_Generator import _Port_Generator
from azapy.Generators.ModelPipeline import ModelPipeline
from azapy.Engines.InvDDEngine import InvDDEngine


class Port_InvDD(_Port_Generator):
    """
    Backtesting Inverse Maximum Drawdown portfolio periodically rebalanced.
    
    **Attributes**
        * `pname` : `str` - portfolio name
        * `ww` : `pandasDataFrame` - portfolio weights at each rebalancing date
        * `port` : `pandas.Series` - portfolio historical time-series
        * `schedule` : `pandas.DataFrame` - rebalancing schedule
       
    The most important method is `set_model`. It must be called before any
    other method.
    """                        
    def set_model(self, hlength=3.25, verbose=False):
        """
        Set model parameters and evaluate the portfolio time-series.
        
        Parameters
        ----------
        hlength : `float`, optional
            The length in year of the historical calibration period relative 
            to `'Dfix'`. A fractional number will be rounded to an integer number 
            of months. The default is `3.25` years. 
        verbose : Boolean, optional
            Sets verbose mode. The default is `False`.

        Returns
        -------
        `pandas.DataFrame` : The portfolio time-series in the format 'date', 
        'pcolname'.
        """
        mod = InvDDEngine(colname=self.col_calib, freq=self.freq,
                          hlength=hlength)
        return super().set_model(ModelPipeline([mod]), verbose)