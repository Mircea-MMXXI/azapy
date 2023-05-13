from .InvVolEngine import InvVolEngine
from azapy.Util.drawdown import max_drawdown

class InvDDEngine(InvVolEngine):
    """
    Inverse maximum drawdown portfolio.
    
    **Attributes**
        * status : `int` - computation status (`0` - success, any other 
          value indicates an error)
        * ww : `pandas.Series` - portfolio weights
        * name : `str` - portfolio name
    """
    def _calc_ww(self):
        self.ww = 1. / self.rrate.apply(lambda x: max_drawdown(x)[0]).abs()
        self.ww /= self.ww.sum(numeric_only=True)
