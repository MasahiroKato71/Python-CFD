from __future__ import annotations

import Initialization 
from PhysicalProperty import *
from Reconstruction import *
from Solver import *
from TimeIntegrationMethod import *
          

lenght = 100
pLow, pHigh = .1, 1.
rhoLow, rhoHigh = .1, 1.
xMin, xMax = -1.0 ,1.0

rho, rhoU, rhoE, p, x = Initialization.ShockTube(lenght, rhoLow, rhoHigh, pLow, pHigh, xMin, xMax, order=2)
riemannRoe = RiemannRoe(epsilon=0.15)
reconstructer = Reconstuctions(Limiter.VanLeer)
rk2 = RungeKutta2(rho, rhoU, rhoE, p, x, riemannRoe, reconstructer.TVD, cflnum=.5)
rk2(starttime=0, endtime=0.55, figName="VanLeer.png")
