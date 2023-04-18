from .InvVolEngine import InvVolEngine

class InvVarEngine(InvVolEngine):
    """
    Inverse variance portfolio.
    
    Methods:
        * getWeights
        * getPositions
        * set_rrate
        * set_mktdata   
    Attributes:
        * status
        * ww
        * name
    """ 
    def _calc_ww(self):
        self.ww = 1. / self.rrate.var(numeric_only=True)
        self.ww /= self.ww.sum(numeric_only=True)