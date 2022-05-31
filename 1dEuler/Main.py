from __future__ import annotations

import Initialization 
from PhysicalProperty import *
from Reconstruction import *
from Solver import *
from TimeIntegrationMethod import *
          

lenght = 100
pLow, pHigh = 0.1, 1.0
rhoLow, rhoHigh = 0.1, 1.0
xMin, xMax = -1.0 ,1.0

rho, rhoU, rhoE, p, x = Initialization.ShockTube(lenght, rhoLow, rhoHigh, pLow, pHigh, xMin, xMax)
riemannRoe = RiemannRoe(rho, rhoU, rhoE, p, epsilon=0.15)
eulearFront = EulearFront(rho, rhoU, rhoE, p, x, riemannRoe, Reconstuctions.piece_wise, cflnum=.5)
eulearFront(starttime=0, endtime=0.55, figName="WindUp.png")
