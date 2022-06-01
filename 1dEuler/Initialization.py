from __future__ import annotations

from PhysicalProperty import *
from Reconstruction import *
from Solver import *


def ShockTube(lenght:int, rhoLow:float, rhoHigh:float, pLow:float, pHigh:float, xMin:float, xMax:float, gamma:float=1.4, order:int=1):
    dx = (xMax - xMin) / lenght
    lenght += 2*order # ゴーストセル分加算
    rho = Rho(lenght)
    rhoU = RhoU(lenght)
    rhoE = RhoE(lenght)
    p = P(gamma)
    x = X(lenght)
    
    x.boundaryValue[0] = xMin - dx
    for i in range(1,lenght+1):
        x.boundaryValue[i] = x.boundaryValue[i-1] + dx
        
    x.SetValue()
        
    for i, x_ in enumerate(x.value):
        if x_ < (xMax+xMin)/2:
            p_ = pHigh
            u_ = 0.
            rho_ = rhoHigh
        else:
            p_ = pLow
            u_ = 0.
            rho_ = rhoLow
        rho.value[i] = rho_
        rhoU.value[i] = rho_ * u_
        rhoE.value[i] = p_/(gamma - 1) - .5 * u_**2 * rho_
        
    return rho, rhoU, rhoE, p, x
