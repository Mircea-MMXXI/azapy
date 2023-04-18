import pandas as pd

from ._Port_Generator import _Port_Generator
from azapy.Generators.ModelPipeline import ModelPipeline
from azapy.Engines.ConstWEngine import ConstWEngine


class Port_ConstW(_Port_Generator):
    """
    Backtesting Constant Weighted portfolio periodically rebalanced.
    
    Methods:
        * set_model
        * get_port
        * get_weights
        * get_nshares
        * get_account
        * get_mktdata
        * port_view
        * port_view_all
        * port_drawdown
        * port_perf
        * port_annual_returns
        * port_monthly_returns
        * port_period_returns
        * port_period_perf
    Attributs:
        * pname
        * ww
        * port
        * schedule
    """                          
    def set_model(self, ww=None, verbose=False):
        """
        Sest model weights.

        Parameters
        ----------
        ww : `pandas.Series`, optional;
            Portfolio weights per symbol. If it is set to `None` then 
            the EWP (Equal Weighted Portfolio) will be considered.
            The default is `None`.
        `verbose` : Boolean, optional;
            Sets verbose mode. The default is `False`.

        Returns
        -------
        `pandas.DataFrame`;
            The portfolio time-series in the format "date", "pcolname".
        """
        if ww is None:
            # set EWP
            return super().set_model(ModelPipeline(["EWP"]), verbose)
        else:
            mod = ConstWEngine(ww)
            return super().set_model(ModelPipeline([mod]), verbose)