import math, matplotlib.pyplot as plt, numpy as np
from scipy.optimize import fsolve
g = 1.4 #gamma

#NORMAL SHOCK RELATIONS
def normal_shock(M1):
    M2 = math.sqrt((M1**2 + 2/(g-1))/(2*g*M1**2/(g-1) - 1))
    P2_P1 = (1+g*M1**2)/(1+g*M2**2)
    R2_R1 = ((g+1)*M1**2)/((g-1)*M1**2+2)
    T2_T1 = (1 + (g-1)*M1**2/2)/(1 + (g-1)*M2**2/2)
    #print(f"Values at M1 = {M1}: \nM2 = {M2}, P2/P1 = {P2_P1}, R2/R1 = {R2_R1}, T2/T1 = {T2_T1}")
    return M2, P2_P1, R2_R1, T2_T1

#OBLIQUE SHOCK RELATIONS
def oblique_shock_angle(theta, M1):
    def equation(x):
        x = x[0]
        return 2/math.tan(x) * (M1**2*(math.sin(x))**2 -1)/(M1**2*(g+math.cos(2*x)) + 2) - math.tan(theta)
    beta, info, ier, msg = fsolve(equation, theta, full_output=True)
    if ier != 1:
        raise RuntimeError("No oblique shock solution")
    return beta[0]
def oblique_shock(theta, M1):
    beta = oblique_shock_angle(theta, M1)
    M2_n, P2_P1, R2_R1, T2_T1 = normal_shock(M1*math.sin(beta))
    M2 = M2_n/(math.sin(beta-theta))
    #print(f"Values at M1 = {M1}: \nbeta = {math.degrees(beta)}, M2_n = {M2_n}, M2 = {M2}, P2/P1 = {P2_P1}, R2/R1 = {R2_R1}, T2/T1 = {T2_T1}")
    return M2, P2_P1, R2_R1, T2_T1

#EXPANSION FAN RELATIONS
def nu(M1):
    return math.sqrt((g+1)/(g-1))*math.atan(math.sqrt((g-1)*(M1**2-1)/(g+1)))-math.atan(math.sqrt(M1**2-1))
def expansion_fan_M2(theta, M1):
    #Find nu(M2)
    nu_M2 = theta + nu(M1)
    #Find M2
    def equation(x):
        x = x[0]
        return math.sqrt((g+1)/(g-1))*math.atan(math.sqrt((g-1)*(x**2-1)/(g+1)))-math.atan(math.sqrt(x**2-1)) - nu_M2
    M2, info, ier, msg = fsolve(equation, 1.5, full_output=True)
    if ier != 1:
        raise RuntimeError("Error in isentropic flow equations")
    return M2[0]
def expansion_fan(theta, M1):
    M2 = expansion_fan_M2(theta, M1)
    T2_T1 = (1 + (g-1)/2*M1**2)/(1 + (g-1)/2*M2**2)
    P2_P1 = (T2_T1)**(g/(g-1))
    R2_R1 = (P2_P1)**(1/g)
    return M2, P2_P1, R2_R1, T2_T1

#DIAMOND SHAPED AIRFOIL
a = 340 #speed of sound in m/s
rho = 1 #density in kg/m^3
p_inf = rho*a**2/g
def airfoil_pres(theta, alpha, M_inf):
    theta = math.radians(theta)
    alpha = math.radians(alpha)
    v_inf = M_inf*a
    
    delta_AD = alpha + theta/2 #always a shock with positive AoA
    M_AD, P2_P1_AD, _, _ = oblique_shock(delta_AD, M_inf)
    p_AD = p_inf * P2_P1_AD
    Cp_AD = (p_AD-p_inf)/(0.5*rho*v_inf**2)

    delta_DC = theta #expansion
    M_DC, P2_P1_DC, _, _ = expansion_fan(abs(delta_DC), M_AD)
    p_DC = p_AD * P2_P1_DC
    Cp_DC = (p_DC-p_inf)/(0.5*rho*v_inf**2)
    
    delta_AB = alpha - theta/2

    if delta_AB >= 0: #expansion or no change
        M_AB, P2_P1_AB, _, _ = expansion_fan(delta_AB, M_inf)
        print(f"expansion, {M_AB, P2_P1_AB}")
    else: #shock
        M_AB, P2_P1_AB, _, _ = oblique_shock(abs(delta_AB), M_inf)
        print("shock")

    p_AB = p_inf * P2_P1_AB
    Cp_AB = (p_AB-p_inf)/(0.5*rho*v_inf**2)

    delta_BC = theta #always an expansion with positive AoA
    M_BC, P2_P1_BC, _, _ = expansion_fan(delta_BC, M_AB)
    p_BC = p_AB * P2_P1_BC
    Cp_BC = (p_BC-p_inf)/(0.5*rho*v_inf**2)
    return Cp_AB, Cp_BC, Cp_AD, Cp_DC

#LIFT AND DRAG COEFFICIENTS
def airfoil_cLcD(af_angle, alpha, M):
    cPs = airfoil_pres(af_angle, alpha, M)
    af_angle = math.radians(af_angle)
    alpha = math.radians(alpha)
    s_c = 1/(2*math.cos(af_angle/2))
    theta = af_angle / 2 # already in radians

    betas = [+theta, math.pi - theta, -theta, math.pi + theta]
    normal_angles = [theta + math.pi/2, (math.pi - theta) - math.pi/2, -theta - math.pi/2, (math.pi + theta) + math.pi/2]
    
    upper = [True, True, False, False]
    cL = 0
    cD = 0
    for normal_angle, cP in zip(normal_angles, cPs):
        cF = -cP * s_c
        dcD = cF * math.cos(normal_angle)
        dcL = cF * math.sin(normal_angle)
        cL += dcL
        cD += dcD
    return cL, cD
def airfoil_cLcD(af_angle_deg, alpha_deg, M):
    cp_AB, cp_BC, cp_AD, cp_DC = airfoil_pres(af_angle_deg, alpha_deg, M)

    theta_deg = af_angle_deg / 2.0
    theta = math.radians(theta_deg)
    alpha = math.radians(alpha_deg)

    deltaAB = theta - alpha
    deltaBC = theta + alpha

    s_c = 1.0 / (2.0 * math.cos(theta))

    CL = s_c * (
    (cp_AD * math.cos(deltaBC) + cp_DC * math.cos(deltaAB)) -
    (cp_AB * math.cos(deltaAB) + cp_BC * math.cos(deltaBC))
    )

    CD = s_c * (
    (cp_AB * math.sin(deltaAB) + cp_AD * math.sin(deltaBC)) -
    (cp_BC * math.sin(deltaBC) + cp_DC * math.sin(deltaAB))
    )

    return CL, CD
alphas = np.linspace(0, 10, 100)
cL = []
cD = []
af_angle = 12
M = 2
cD0 = airfoil_cLcD(af_angle, 0, M)[1]
for alpha in alphas:
    cL.append(airfoil_cLcD(af_angle, alpha, M)[0])
    cD.append(airfoil_cLcD(af_angle, alpha, M)[1])
plt.plot(alphas, cL)
plt.title("cL vs alpha at M=2")
plt.show()
plt.plot(alphas, cD)
plt.title("cD vs alpha at M=2")
plt.show()