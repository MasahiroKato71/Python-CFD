from __future__ import annotations

import copy

from ModelBase import *
from PhysicalProperty import *
from Reconstruction import *
from Solver import *

    
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
         