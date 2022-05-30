from __future__ import annotations

from typing import Final

import matplotlib.pyplot as plt

from PhysicalProperty import *
from Reconstruction import *
from Solver import *


def update(rho:Rho, rhoU:RhoU, rhoE:RhoE, x:X):
     if not (type(rho) is Rho and type(rhoU) is RhoU and type(rhoE) is RhoE):
          raise TypeError
     if type(x) is not X:
          raise TypeError
     for i in range(1, len(rho.value)):
          rho.value[i] += dt*(rho.f[i] - rho.f[i+1]) / (x.value[i+1] - x.value[i])
          rhoU.value[i] += dt*(rhoU.f[i] - rhoU.f[i+1]) / (x.value[i+1] - x.value[i])
          rhoE.value[i] += dt*(rhoE.f[i] - rhoE.f[i+1]) / (x.value[i+1] - x.value[i])
          
          
def cfl(rho:Rho, rhoU:RhoU, rhoE:RhoE, p:P, x:X, cflnum:float=0.5) -> float:
     if not (type(rho) is Rho and type(rhoU) is RhoU and type(rhoE) is RhoE):
          raise TypeError
     if not (type(p) is P and type(x) is X and type(cflnum) is float):
          raise TypeError
     dt = 1e2
     for i in range(len(rho.value)):
          u = rhoU.value[i] / rho.value[i]
          c = math.sqrt(p.gamma * p(rho.value[i], rhoU.value[i], rhoE.value[i])/rho.value[i])
          lambda1 = abs(u)
          lambda2 = abs(u + c)
          lambda3 = abs(u - c)
          lambdaMax = max(lambda1, lambda2, lambda3, 0.1)
          dt = min(dt, cflnum*(x.value[i+1] - x.value[i])/lambdaMax)
     
     return dt
          
          
LENGTH:Final[int] = 100
XMIN:Final[float] = -1.
XMAX:Final[float] = 1.
TSTOP:Final[float] = .5

GAMMA:Final[float] = 1.4
CFL:Final[float] = .5

pltInterval = 25

rho = Rho(LENGTH)
rhoU = RhoU(LENGTH)
rhoE = RhoE(LENGTH)
p = P(GAMMA)

x = X(LENGTH)
dx = (XMAX - XMIN)/LENGTH

# 初期化処理
x.value[0] = XMIN
for i in range(1,len(x.value)):
     x.value[i] = x.value[i-1] + dx
     
for i, x_ in enumerate(x.value[:-1]):
     if x_ < 0:
          p_ = 1.
          u_ = 0.
          rho_ = 1.
     else:
          p_ = 0.1
          u_ = 0.
          rho_ = 0.1
     rho.value[i] = rho_
     rhoU.value[i] = rho_ * u_
     rhoE.value[i] = p_/(GAMMA - 1) - .5 * u_**2 * rho_
               

t = 0.
n = 0.

riemannRoe = RiemannRoe(rho, rhoU, rhoE, p, epsilon=0.15)
plt.figure(figsize=(16,9))
plt.rcParams["font.size"] = 26
plt.plot(x.value[:-1], rho.value, lw=2, label="t=0")

while t < TSTOP:
     n += 1
     dt = cfl(rho, rhoU, rhoE, p, x, CFL)
     t += dt
     Reconstuctions.piece_wise(rho)
     Reconstuctions.piece_wise(rhoU)
     Reconstuctions.piece_wise(rhoE)
     riemannRoe()
     update(rho, rhoU, rhoE, x)
     if n % pltInterval == 0:
          plt.plot(x.value[:-1], rho.value, lw=2, label=f"t={t:.2f}")

plt.xlabel("x")
plt.ylabel("rho")
plt.legend()
plt.savefig("Eulear.png")
     
     
     