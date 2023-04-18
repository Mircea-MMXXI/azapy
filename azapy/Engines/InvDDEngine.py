from .InvVolEngine import InvVolEngine
from azapy.util.drawdown import max_drawdown

class InvDDEngine(InvVolEngine):
    """
    Inverse maximum drwadown portfolio.
    
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
        self.ww = 1. / self.rrate.apply(lambda x: max_drawdown(x)[0]).abs()
        self.ww /= self.ww.sum(numeric_only=True)
