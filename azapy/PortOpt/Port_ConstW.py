from ._Port_Generator import _Port_Generator
from azapy.Generators.ModelPipeline import ModelPipeline
from azapy.Engines.ConstWEngine import ConstWEngine


class Port_ConstW(_Port_Generator):
    """
    Backtesting Constant Weighted portfolio periodically rebalanced.
    
    **Attributes**
        * `pname` : `str` - portfolio name
        * `ww` : `pandasDataFrame` - portfolio weights at each rebalancing date
        * `port` : `pandas.Series` - portfolio historical time-series
        * `schedule` : `pandas.DataFrame` - rebalancing schedule
       
    The most important method is `set_model`. It must be called before any
    other method.
    """                            
    def set_model(self, ww=None, verbose=False):
        """
        Sets model weights.

        Parameters
        ----------
        ww : `pandas.Series`, optional
            Portfolio weights per symbol. If it is set to `None` then 
            the EWP (Equal Weighted Portfolio) will be considered.
            The default is `None`.
        verbose : Boolean, optional
            Sets verbose mode. The default is `False`.

        Returns
        -------
        `pandas.DataFrame` : Portfolio time-series in the format 'date', 
        'pcolname'.
        """
        if ww is None:
            # set EWP
            return super().set_model(ModelPipeline(["EWP"]), verbose)
        else:
            mod = ConstWEngine(ww)
            return super().set_model(ModelPipeline([mod]), verbose)