from __future__ import annotations

import copy

import matplotlib.pyplot as plt

from ModelBase import *
from PhysicalProperty import *
from Reconstruction import *
from Solver import *


# ToDo:別の場所の親クラスを呼び出したときにRhoなどでNameErrorが発生するため移設中
class TimeIntegrationMethodBase():
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
        self.t = None
        
    def __call__(self, starttime, endtime, figName=None, pltInterval=25):
        self.t = starttime
        n = 0

        if figName:
            self.SetupFig()

        while self.t < endtime:
            n += 1
            
            dt = self.dtCalc()
            
            self.Update(dt)
            
            self.t += dt
            
            if n % pltInterval == 0:
                self.Plot()

        if figName:
            self.OutputFig(figName)
            
    def Update(self, dt:float):
        raise NotImplementedError
    
    def dtCalc(self) -> float:
        dt = 1e2
        for i in range(1,len(self.x.value)-1):
            u = self.rhoU.value[i] / self.rho.value[i]
            c = math.sqrt(self.p.gamma * self.p(self.rho.value[i], self.rhoU.value[i], self.rhoE.value[i])/ self.rho.value[i])
            lambda1 = abs(u)
            lambda2 = abs(u + c)
            lambda3 = abs(u - c)
            lambdaMax = max(lambda1, lambda2, lambda3, 0.1)
            dt = min(dt, self.cflnum*(self.x.boundaryValue[i+1] - self.x.boundaryValue[i])/lambdaMax)
        
        return dt
    
    def SetupFig(self) -> None:
        plt.rcParams["font.size"] = 18
        self.fig, (self.figRho, self.figP) = plt.subplots(figsize=(16,12), nrows=2, ncols=1, sharex=True)
        self.Plot()

    def Plot(self) -> None:
        self.figRho.plot(self.x.value[self.solver.order:-1*self.solver.order], self.rho.value[self.solver.order:-1*self.solver.order], label=f"t={self.t:.2f}")
        self.figP.plot(self.x.value[self.solver.order:-1*self.solver.order],\
            [self.p(self.rho.value[i], self.rhoU.value[i], self.rhoE.value[i]) for i in range(self.solver.order, len(self.x.value) - self.solver.order)], label=f"t={self.t:.2f}")

    def OutputFig(self, figName) -> None:
        self.Plot()
        self.figRho.set_ylabel("rho")
        self.figRho.legend()
        self.figP.set_xlabel("x")
        self.figP.set_ylabel("p")
        self.figP.legend()
        plt.savefig(figName)
    
    
class EulearFront(TimeIntegrationMethodBase):
    def __init__(self, rho:Rho, rhoU:RhoU, rhoE:RhoE, p:P, x:X, solver:SolverBase, reconstructer:function, cflnum:float):
        super().__init__(rho, rhoU, rhoE, p, x, solver, reconstructer, cflnum)
        self.solver.order = 1
        
    def Update(self, dt:float) -> None:
        if type(dt) is not float:
            raise TypeError
        
        self.reconstructer(self.rho)
        self.reconstructer(self.rhoU)
        self.reconstructer(self.rhoE)
        self.solver()
        for i in range(1, len(self.x.value)):
            self.rho.value[i] += dt*(self.rho.f[i] - self.rho.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i])
            self.rhoU.value[i] += dt*(self.rhoU.f[i] - self.rhoU.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i])
            self.rhoE.value[i] += dt*(self.rhoE.f[i] - self.rhoE.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i])
            
            
class RungeKutta2(TimeIntegrationMethodBase):
    def __init__(self, rho:Rho, rhoU:RhoU, rhoE:RhoE, p:P, x:X, solver:SolverBase, reconstructer:function, cflnum:float):
        super().__init__(rho, rhoU, rhoE, p, x, solver, reconstructer, cflnum)
        self.solver.order = 2
            
    def Update(self, dt:float) -> None:
        if type(dt) is not float:
            raise TypeError
        
        self.reconstructer(self.rho)
        self.reconstructer(self.rhoU)
        self.reconstructer(self.rhoE)
        self.solver()
        
        rhoAster, rhoUAster, rhoEAster = copy.copy(self.rho.value), copy.copy(self.rhoU.value), copy.copy(self.rhoE.value)
        LhRho, LhRhoU, LhRhoE = [0. for _ in range(len(self.x.value))], [0. for _ in range(len(self.x.value))], [0. for _ in range(len(self.x.value))]
        for i in range(1, len(self.x.value)):
            LhRho[i] = (self.rho.f[i] - self.rho.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i])
            LhRhoU[i] = (self.rhoU.f[i] - self.rhoU.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i])
            LhRhoE[i] = (self.rhoE.f[i] - self.rhoE.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i])
            rhoAster[i] += dt*LhRho[i]
            rhoUAster[i] += dt*LhRhoU[i]
            rhoEAster[i] += dt*LhRhoE[i]
            
        rhoValue, rhoUValue, rhoEValue = copy.copy(self.rho.value), copy.copy(self.rhoU.value), copy.copy(self.rhoE.value)
        self.rho.value, self.rhoU.value, self.rhoE.value = rhoAster, rhoUAster, rhoEAster
                
        self.reconstructer(self.rho)
        self.reconstructer(self.rhoU)
        self.reconstructer(self.rhoE)
        self.solver()
        
        for i in range(len(self.x.value)):
            self.rho.value[i] = rhoValue[i] + .5 * dt*(LhRho[i] + (self.rho.f[i] - self.rho.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i]))
            self.rhoU.value[i] = rhoUValue[i] + .5 * dt*(LhRhoU[i] + (self.rhoU.f[i] - self.rhoU.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i]))
            self.rhoE.value[i] = rhoEValue[i] + .5 * dt*(LhRhoE[i] + (self.rhoE.f[i] - self.rhoE.f[i+1]) / (self.x.boundaryValue[i+1] - self.x.boundaryValue[i]))
         