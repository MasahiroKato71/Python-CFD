from __future__ import annotations

import math

import matplotlib.pyplot as plt

import PhysicalProperty


class PhysicalPropertyBase():
    def __init__(self, length:int):
        if type(length) is not int:
            raise TypeError
        
        self.value = [0. for _ in range(length)]
        self.leftBoundaryValue = [0. for _ in range(length+1)]
        self.rightBoundaryValue = [0. for _ in range(length+1)]
        self.f = [0. for _ in range(length+1)]
        
    def F(self) -> float:
        raise NotImplementedError
    
    
class SolverBase():
    def __call__() -> None:
        raise NotImplementedError


# TimeIntegrationMethodBase:循環参照が発生するためTimeIntegratuonMethod.pyに移設中
class TimeIntegrationMethodBase():
    def __init__(self, rho:PhysicalProperty.Rho, rhoU:PhysicalProperty.RhoU, rhoE:PhysicalProperty.RhoE, p:PhysicalProperty.P, x:PhysicalProperty.X, solver:SolverBase, reconstructer:function, cflnum:float):
        if not (type(rho) is PhysicalProperty.Rho and type(rhoU) is PhysicalProperty.RhoU and type(rhoE) is PhysicalProperty.RhoE \
            and type(p) is PhysicalProperty.P and type(x) is PhysicalProperty.X):
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
