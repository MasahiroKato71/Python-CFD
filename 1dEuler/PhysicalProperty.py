from __future__ import annotations

import ModelBase
from typing import Final


# セル位置モデル
class X():
    def __init__(self, length:int):
        if type(length) is not int:
            raise TypeError
        
        self.value = [0. for _ in range(length)]
        self.boundaryValue = [0. for _ in range(length+1)]
        
    def SetValue(self):
        for i in range(len(self.value)-1):
            self.value[i] = (self.boundaryValue[i] + self.boundaryValue[i+1]) / 2

# 圧力モデル
class P():
    def __init__(self, gamma:float=1.4):
        if type(gamma) is not float:
            raise TypeError
        
        self.gamma:Final[float] = gamma
        
    def __call__(self, rho:float, rhoU:float, rhoE:float) -> float:
        if not (type(rho) is float and type(rhoU) is float and type(rhoE) is float):
            raise TypeError
        
        return (self.gamma - 1) * (rhoE - rhoU**2 / rho)
    
    
# 密度モデル（質量保存則）
class Rho(ModelBase.PhysicalPropertyBase):
    @staticmethod
    def F(rhoU) -> float:
        if type(rhoU) is not float:
            raise TypeError
        
        return rhoU


# 運動量モデル（運動量保存則）
class RhoU(ModelBase.PhysicalPropertyBase):
    @staticmethod
    def F(rho:float, rhoU:float, p:float) -> float:
        if not (type(rho) is float and type(rhoU) is float and type(p) is float):
            raise TypeError
        
        return rhoU ** 2 / rho + p


# 総エネルギーモデル（エネルギー保存則）
class RhoE(ModelBase.PhysicalPropertyBase):
    @staticmethod
    def F(rho:float, rhoU:float, rhoE:float, p:float) -> float:
        if not (type(rho) is float and type(rhoU) is float and type(rhoE) is float and type(p) is float):
            raise TypeError
        
        return (rhoE + p) * rhoU / rho
    