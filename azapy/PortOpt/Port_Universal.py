import pandas as pd
from azapy.Engines.UniversalEngine import UniversalEngine
from azapy.Generators.Port_Generator import Port_Generator

class Port_Universal(Port_Generator):
    """
    Backtesting Cover's (1996) Universal portfolio.
    
    **Attributes**
        * `pname` : `str` - portfolio name
        * `ww` : `pandasDataFrame` - portfolio weights at each rebalancing date
        * `port` : `pandas.Series` - portfolio historical time-series
        * `schedule` : `pandas.DataFrame` - rebalancing schedule
       
    The most important method is `set_model`. It must be called before any
    other method.
    
    Note: if '_CASH_' asset is not explicitly present in the `mktdata` (as
    a portfolio cash component), then it will be added to the portfolio 
    with a position (weight) net `0`.
    """                  
    def set_model(self, mc_paths=100, nr_batches=16, 
                  variance_reduction=True, dirichlet_alpha=None, 
                  mc_seed=None, verbose=False):
        """
        Sets model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        mc_paths : positive `int`, optional
            Number of Monte Carlo simulations
            The default is `100`.
        nr_batches : positive `int`, optional
            Number of Monte Carlo simulation batches. Each batch runs
            `mc_paths` simulations. The computation is multithreaded for 
            `nr_batches > 1`.
            The default is `16`.
        variance_reduction : Boolean, optional
            If set to `True`, then the antithetic variance reduction based on 
            all possible permutations of the basket components is deployed.
            In this case the total number of MC simulations is
            `mc_paths * nr_batches * factorial(number of portfolio components)`.
            For example, a portfolio of 6 assets, with default values for 
            `mc_paths` and `nr_batches` implies `100 * 16 * 720 = 1,1520,000`
            MC simulations. If it is set to `False`, then the total number
            of simulations is `mc_paths * nr_batches`.
            The default is `True`.
        mc_seed : positive `int`, optional
            Random number generator seed. If it set to a positive `int` 
            the Monte Carlo random number generator is reseeded.
            The default is `None`.
        dirichlet_alpha : `list`, optional
            The alpha values for Dirichlet random vector generator. Must have
            the size equal to the number of portfolio components. If it is
            set to `None` then the uniform random generator vectors 
            in a n-simplex is used (equivalent to Flat Dirichlet random 
            vector generator, i.e., alpha=[1] * n).
            The default is `None`.
        verbose : Boolean, optional
            Logical flag triggering the verbose mode. The default is `False`.

        Returns
        -------
        `pandas.DataFrame` : The portfolio time-series in the format 'date', 
        'pcolname'.
        """
        self._set_schedule()
        self.wwModel = UniversalEngine(mktdata=self.mktdata, 
                                       schedule=self.schedule,
                                       dirichlet_alpha=dirichlet_alpha)
        ww = self.wwModel.getWeights(mc_paths=mc_paths,
                                     nr_batches=nr_batches,
                                     variance_reduction=variance_reduction,
                                     mc_seed=mc_seed,
                                     verbose=verbose)
        
        self.ww = pd.merge(self.schedule, ww, left_on='Dfix', right_index=True)
        if '_CASH_' not in self.ww.columns:
            self.ww['_CASH_'] = 0
        self._port_calc()

        return self.port
