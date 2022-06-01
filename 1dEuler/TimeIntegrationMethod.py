from __future__ import annotations

from typing import Final

import matplotlib.pyplot as plt

from PhysicalProperty import *
from Reconstruction import *
from Solver import *


class EulearFront(TimeIntegrationMethodBase):
    def __init__(self, rho:Rho, rhoU:RhoU, rhoE:RhoE, p:P, x:X, solver:SolverBase, reconstructer:function, cflnum:float):
        if not (type(rho) is Rho and type(rhoU) is RhoU and type(rhoE) is RhoE and type(p) is P and type(x) is X):
            raise TypeError
        if not (isinstance(solver, SolverBase) and callable(reconstructer)):
            raise TypeError
        if type(cflnum) is not float:
            raise TypeError
                 
        self.rho = rho
        self.rhoU = rhoU
        self.rhoE = rhoE
        self.p = p
        self.x = x
        self.solver = solver
        self.solver.rho = rho
        self.solver.rhoU = rhoU
        self.solver.rhoE = rhoE
        self.solver.p = p
        self.reconstructer = reconstructer        
        self.cflnum = cflnum
    
    def __call__(self, starttime, endtime, figName=None, pltInterval=25):
        t = starttime
        n = 0

        if figName:
            plt.rcParams["font.size"] = 26
            fig, (figRho, figP) = plt.subplots(figsize=(16,12), nrows=2, ncols=1, sharex=True)
            figRho.plot(self.x.value[1:-1], self.rho.value[1:-1], label="t=0.00")
            figP.plot(self.x.value[1:-1], [self.p(self.rho.value[i], self.rhoU.value[i], self.rhoE.value[i]) for i in range(1, len(self.x.value)-1)], label="t=0.00")

        while t < endtime:
            n += 1
            
            dt = self.dtCalc()
            
            self.reconstructer(self.rho)
            self.reconstructer(self.rhoU)
            self.reconstructer(self.rhoE)
            self.solver()
            self.Update(dt)
            
            t += dt
            
            if n % pltInterval == 0:
                figRho.plot(self.x.value[1:-1], self.rho.value[1:-1], label=f"t={t:.2f}")
                figP.plot(self.x.value[1:-1], [self.p(self.rho.value[i], self.rhoU.value[i], self.rhoE.value[i]) for i in range(1, len(self.x.value)-1)], label=f"t={t:.2f}")

        if figName:
            figRho.plot(self.x.value[1:-1], self.rho.value[1:-1], label=f"t={t:.2f}")
            figRho.set_ylabel("rho")
            figRho.legend()
            figP.plot(self.x.value[1:-1], [self.p(self.rho.value[i], self.rhoU.value[i], self.rhoE.value[i]) for i in range(1, len(self.x.value)-1)], label=f"t={t:.2f}")
            figP.set_xlabel("x")
            figP.set_ylabel("p")
            figP.legend()
            plt.savefig(figName)
        
    def Update(self, dt:float) -> None:
        if type(dt) is not float:
            raise TypeError
        for i in range(1, len(self.rho.value)):
            self.rho.value[i] += dt*(self.rho.f[i] - self.rho.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i])
            self.rhoU.value[i] += dt*(self.rhoU.f[i] - self.rhoU.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i])
            self.rhoE.value[i] += dt*(self.rhoE.f[i] - self.rhoE.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i])
          
        
    def dtCalc(self) -> float:
        dt = 1e2
        for i in range(len(self.rho.value)):
            u = self.rhoU.value[i] / self.rho.value[i]
            c = math.sqrt(self.p.gamma * self.p(self.rho.value[i], self.rhoU.value[i], self.rhoE.value[i])/ self.rho.value[i])
            lambda1 = abs(u)
            lambda2 = abs(u + c)
            lambda3 = abs(u - c)
            lambdaMax = max(lambda1, lambda2, lambda3, 0.1)
            dt = min(dt, self.cflnum*(self.x.boundaryValue[i+1] - self.x.boundaryValue[i])/lambdaMax)
        
        return dt
    