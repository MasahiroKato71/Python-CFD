from __future__ import annotations

import math

from ModelBase import SolverBase
from PhysicalProperty import *


class RiemannRoe(SolverBase):
    def __init__(self, epsilon:float=0.15):
        if type(epsilon) is not float:
            raise TypeError
        
        self.rho = self.rhoU = self.rhoE = self.p = None
        self.order = None
        self.epsilon = epsilon
        
    def __call__(self) -> None:
        for i in range(self.order, len(self.rho.leftBoundaryValue)-self.order):
            # プリミティブ変数の算出
            rhoL = self.rho.leftBoundaryValue[i]
            rhoR = self.rho.rightBoundaryValue[i]
            uL = self.rhoU.leftBoundaryValue[i] / self.rho.leftBoundaryValue[i]
            uR = self.rhoU.rightBoundaryValue[i] / self.rho.rightBoundaryValue[i]
            pL = self.p(self.rho.leftBoundaryValue[i], self.rhoU.leftBoundaryValue[i], self.rhoE.leftBoundaryValue[i])
            pR = self.p(self.rho.rightBoundaryValue[i], self.rhoU.rightBoundaryValue[i], self.rhoE.rightBoundaryValue[i])
            hL = (self.rhoE.leftBoundaryValue[i] + pL) / self.rho.leftBoundaryValue[i]
            hR = (self.rhoE.rightBoundaryValue[i] + pR) / self.rho.rightBoundaryValue[i]
            eL = self.rhoE.leftBoundaryValue[i] / self.rho.leftBoundaryValue[i]
            eR = self.rhoE.rightBoundaryValue[i] / self.rho.rightBoundaryValue[i]
                        
            # Roe平均の算出
            rhoAve = math.sqrt(rhoL * rhoR)
            uAve = (uL * math.sqrt(rhoL) + uR * math.sqrt(rhoR)) / (math.sqrt(rhoL) + math.sqrt(rhoR))
            hAve = (hL * math.sqrt(rhoL) + hR * math.sqrt(rhoR)) / (math.sqrt(rhoL) + math.sqrt(rhoR))
            cAve = math.sqrt((self.p.gamma - 1)*(hAve - uAve**2/2))
            
            # 特性速度の算出
            lambda1 = Harten(uAve, self.epsilon)
            lambda2 = Harten(uAve + cAve, self.epsilon)
            lambda3 = Harten(uAve - cAve, self.epsilon)
            
            # 数値流束の算出
            dw1 = rhoR - rhoL - (pR - pL)/cAve**2
            dw2 = uR - uL + (pR - pL) / (rhoAve*cAve)
            dw3 = uR - uL - (pR - pL) / (rhoAve*cAve)
            
            self.rho.f[i] = .5 * (self.rho.F(rhoR*uR) + self.rho.F(rhoL*uL)) \
                - .5 * (lambda1*dw1 + lambda2*rhoAve/(2*cAve)*dw2 - lambda3*rhoAve/(2*cAve)*dw3)
            self.rhoU.f[i] = .5 * (self.rhoU.F(rhoR, rhoR*uR, pR) + self.rhoU.F(rhoL, rhoL*uL, pL)) \
                - .5 * (lambda1*dw1*uAve + lambda2*rhoAve/(2*cAve)*dw2*(uAve+cAve) - lambda3*rhoAve/(2*cAve)*dw3*(uAve-cAve))
            self.rhoE.f[i] = .5 * (self.rhoE.F(rhoR, rhoR*uR, rhoR*eR, pR) + self.rhoE.F(rhoL, rhoL*uL, rhoL*eL, pL)) \
                - .5 * (lambda1*dw1*(uAve**2/2) + lambda2*rhoAve/(2*cAve)*dw2*(hAve+cAve*uAve) - lambda3*rhoAve/(2*cAve)*dw3*(hAve-cAve*uAve))

        # 境界条件
        for i in range(self.order):
            self.rho.f[i] = self.rho.f[self.order]
            self.rhoU.f[i] = self.rhoU.f[self.order]
            self.rhoE.f[i] = self.rhoE.f[self.order]
        
            self.rho.f[-1*i - 1] = self.rho.f[-1*self.order - 1]
            self.rhoU.f[-1*i - 1] = self.rhoU.f[-1*self.order - 1]
            self.rhoE.f[-1*i - 1] = self.rhoE.f[-1*self.order - 1]
        

def Harten(alpha:float, epsilon:float) -> float:
    if not (type(alpha) is float and type(epsilon) is float):
        raise TypeError
    
    if abs(alpha) < 2*epsilon:
        return alpha**2 / (4 * epsilon) + epsilon       
    else:
        return abs(alpha)
    