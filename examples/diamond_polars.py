import numpy as np
import matplotlib.pyplot as plt
from compressible_flow.airfoil import airfoil_cLcD

af_angle = 12
M = 2
alphas = np.linspace(0, 10, 100)
cL = []
cD = []

for alpha in alphas:
    cL.append(airfoil_cLcD(af_angle, alpha, M)[0])
    cD.append(airfoil_cLcD(af_angle, alpha, M)[1])
plt.plot(alphas, cL)
plt.title("CL vs alpha at M=2")
plt.xlabel("alpha (deg)")
plt.ylabel("CL")
plt.show()

plt.plot(alphas, cD)
plt.title("CD vs alpha at M=2")
plt.xlabel("alpha (deg)")
plt.ylabel("CD")
plt.show()