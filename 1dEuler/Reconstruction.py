from __future__ import annotations
from ModelBase import *


class Reconstuctions():
    @staticmethod
    def piece_wise(physicalProperty:PhysicalPropertyBase):
        for i in range(1,len(physicalProperty.value)):
            physicalProperty.leftBoundaryValue[i] = physicalProperty.value[i-1]
            physicalProperty.rightBoundaryValue[i] = physicalProperty.value[i]