import os

__all__ = [_[:-3] for _ in os.listdir("src/subs") if _.endswith(".py") and _ != "__init__.py"]

from subs.CN import CN
from subs.SCN import SCN
from subs.NNH import NNH
from subs.NNN import NNN
from subs.NNOH import NNOH
from subs.SSCH3 import SSCH3
from subs.CCCH3 import CCCH3

