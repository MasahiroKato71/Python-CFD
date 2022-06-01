from __future__ import annotations

from ModelBase import *
import Utils

class Reconstuctions():
    def __init__(self, limiter:function = None):
        if not (callable(limiter) or limiter is None):
            raise TypeError
        
        self.limiter = limiter
        
    @staticmethod
    def PieceWise(physicalProperty:PhysicalPropertyBase) -> None:
        if not isinstance(physicalProperty, PhysicalPropertyBase):
            raise TypeError
        
        for i in range(1,len(physicalProperty.value)):
            physicalProperty.leftBoundaryValue[i] = physicalProperty.value[i-1]
            physicalProperty.rightBoundaryValue[i] = physicalProperty.value[i]
            
    def TVD(self, physicalProperty:PhysicalPropertyBase) -> None:
        if not isinstance(physicalProperty, PhysicalPropertyBase):
            raise TypeError
        if not callable(self.limiter):
            raise TypeError
        
        for i in range(2, len(physicalProperty.leftBoundaryValue)-2):
            physicalProperty.leftBoundaryValue[i] = physicalProperty.value[i-1] \
                + .5 * self.limiter((physicalProperty.value[i-1] - physicalProperty.value[i-2]) / Utils.MaxAbs((physicalProperty.value[i] - physicalProperty.value[i-1]), 1e-6)) \
                    * (physicalProperty.value[i] - physicalProperty.value[i-1])
            
            physicalProperty.rightBoundaryValue[i] = physicalProperty.value[i] \
                - .5 * self.limiter((physicalProperty.value[i+1] - physicalProperty.value[i]) / Utils.MaxAbs((physicalProperty.value[i] - physicalProperty.value[i-1]), 1e-6)) \
                    * (physicalProperty.value[i] - physicalProperty.value[i-1])


class Limiter():
    @staticmethod
    def Minmod(r:float) -> float:
        if type(r) is not float:
            raise TypeError
        
        return max(0, min(1, r))
    
    @staticmethod
    def Superbee(r:float) -> float:    
        if type(r) is not float:
            raise TypeError
        
        return max(0, min(1, 2*r), min(2, r))
    
    @staticmethod
    def VanLeer(r:float) -> float:
        if type(r) is not float:
            raise TypeError
        
        return (r + abs(r)) / (1 + abs(r))
    
    @staticmethod
    def VanAlbada(r:float) -> float:
        if type(r) is not float:
            raise TypeError
        
        return (r + r**2) / (1 + r**2)