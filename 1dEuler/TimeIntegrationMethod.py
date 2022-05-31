from __future__ import annotations

from typing import Final

import matplotlib.pyplot as plt

from PhysicalProperty import *
from Reconstruction import *
from Solver import *


class EulearFront():
    def __init__(self, rho:Rho, rhoU:RhoU, rhoE:RhoE, p:P, x:X, solver:SolverBase, reconstructer:function, cflnum:float):
        self.rho = rho
        self.rhoU = rhoU
        self.rhoE = rhoE
        self.p = p
        self.x = x
        self.solver = solver
        self.reconstructer = reconstructer        
        self.cflnum = cflnum
    
    def __call__(self, starttime, endtime, figName=None, pltInterval=25):
        t = starttime
        n = 0

        if figName:
            plt.figure(figsize=(16,9))
            plt.rcParams["font.size"] = 26
            plt.plot(self.x.value[1:-2], self.rho.value[1:-1], lw=2, label="t=0")

        while t < endtime:
            n += 1
            
            dt = self.dtCalc()
            
            self.reconstructer(self.rho)
            self.reconstructer(self.rhoU)
            self.reconstructer(self.rhoE)
            self.solver()
            self.update(dt)
            
            t += dt
            
            if n % pltInterval == 0:
                plt.plot(self.x.value[1:-2], self.rho.value[1:-1], lw=2, label=f"t={t:.2f}")

        if figName:
            plt.plot(self.x.value[1:-2], self.rho.value[1:-1], lw=2, label=f"t={t:.2f}")
            plt.xlabel("x")
            plt.ylabel("rho")
            plt.legend()
            plt.savefig(figName)
        
    def update(self, dt:float) -> None:
        if type(dt) is not float:
            raise TypeError
        for i in range(1, len(self.rho.value)):
            self.rho.value[i] += dt*(self.rho.f[i] - self.rho.f[i+1]) / (self.x.value[i+1] - self.x.value[i])
            self.rhoU.value[i] += dt*(self.rhoU.f[i] - self.rhoU.f[i+1]) / (self.x.value[i+1] - self.x.value[i])
            self.rhoE.value[i] += dt*(self.rhoE.f[i] - self.rhoE.f[i+1]) / (self.x.value[i+1] - self.x.value[i])
          
          
    def dtCalc(self) -> float:
        dt = 1e2
        for i in range(len(self.rho.value)):
            u = self.rhoU.value[i] / self.rho.value[i]
            c = math.sqrt(self.p.gamma * self.p(self.rho.value[i], self.rhoU.value[i], self.rhoE.value[i])/ self.rho.value[i])
            lambda1 = abs(u)
            lambda2 = abs(u + c)
            lambda3 = abs(u - c)
            lambdaMax = max(lambda1, lambda2, lambda3, 0.1)
            dt = min(dt, self.cflnum*(self.x.value[i+1] - self.x.value[i])/lambdaMax)
        
        return dt




     
     
     