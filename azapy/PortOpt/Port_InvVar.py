from ._Port_Generator import _Port_Generator
from azapy.Generators.ModelPipeline import ModelPipeline
from azapy.Engines.InvVarEngine import InvVarEngine


class Port_InvVar(_Port_Generator):
    """
    Backtesting Invese Variance portfolio periodically rebalanced.
    
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
    def set_model(self, hlength=3.25, verbose=False):
        """
        Set model parameters and evaluate the portfolio time-series.
        
        Parameters
        ----------
        `hlength` : `float`, optional;
            The length in year of the historical calibration period relative 
            to `'Dfix'`. A fractional number will be rounded to an integer number 
            of months. The default is `3.25` years. 
        `verbose` : Boolean, optional;
            Sets verbose mode. The default is `False`.

        Returns
        -------
        `pandas.DataFrame`;
            The portfolio time-series in the format "date", "pcolname".
        """
        mod = InvVarEngine(colname=self.col_calib, freq=self.freq,
                           hlength=hlength)
        return super().set_model(ModelPipeline([mod]), verbose)