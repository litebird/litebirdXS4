import numpy as np
from .config import NOISEDIR



class NoiseModel:
    
    
    def __init__(self,nside=512):
        pass
    
    def noise_freq(self,freq,idx):
        raise NotImplementedError