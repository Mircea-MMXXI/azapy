__version__ = "1.2.2"
from .MkT import *
from .Analyzers import *
from .Engines import *
from .Selectors import *
from .PortOpt import *
from .Generators import *
from .Util import *

def version():
    """Returns **azapy** package version"""
    return __version__
