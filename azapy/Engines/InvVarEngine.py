from .InvVolEngine import InvVolEngine

class InvVarEngine(InvVolEngine):
    """
    Inverse variance portfolio.
    
    **Attributes**
        * status : `int` - computation status (`0` - success, any other 
          value indicates an error)
        * ww : `pandas.Series` - portfolio weights
        * name : `str` - portfolio name
    """
    def _calc_ww(self):
        self.ww = 1. / self.rrate.var(numeric_only=True)
        self.ww /= self.ww.sum(numeric_only=True)