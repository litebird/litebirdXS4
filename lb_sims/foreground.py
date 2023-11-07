import healpy as hp
import numpy as np
from .config import FGDIR


class Foregrounds:


    def __init__(self,nside,complexity='low'):
        self.nside = nside

        if complexity not in ['low','medium','high']:
            raise ValueError 'wrong complexity, try low,medium,high!'
        self.models = self.__select_FG__.(complexity)



    def __selectFG__(self,which):
        return {
                'low': ['d9','s4','f1','a1','co1'],
                'medium' : ['d10','s5','f1','a1','co3'],
                'high': ['d12','s7','f1','a2','co3']
                }[which]

    def TQU(self,freq):
        raise NotImplementedError


class CompSep:

    def __init__(self,method='HILC'):
        pass
