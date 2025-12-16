import math
from scipy.optimize import fsolve
from .constants import DEFAULT_GAMMA

def normal_shock(M1: float, g: float = DEFAULT_GAMMA):
    """
    Normal shock relations
    
    :param M1: upstream Mach number
    :type M1: float
    :param g: gamma
    :type g: float
    :return: downstream Mach number, pressure ratio, density ratio, temperature ratio
    :rtype: tuple[float, float, float, float]
    """
    M2 = math.sqrt((M1**2 + 2/(g-1))/(2*g*M1**2/(g-1) - 1))
    P2_P1 = (1+g*M1**2)/(1+g*M2**2)
    R2_R1 = ((g+1)*M1**2)/((g-1)*M1**2+2)
    T2_T1 = (1 + (g-1)*M1**2/2)/(1 + (g-1)*M2**2/2)
    return M2, P2_P1, R2_R1, T2_T1

def oblique_shock_angle(theta: float, M1: float, g: float = DEFAULT_GAMMA) -> float:
    """
    Solves theta-beta-M relation for beta in radians
    
    :param theta: flow deflection angle in [radians]
    :type theta: float
    :param M1: upstream Mach number
    :type M1: float
    :param g: gamma
    :type g: float
    :return: beta - shock wave angle [radians]
    :rtype: float
    """
    def equation(x):
        beta = x[0]
        return 2/math.tan(beta) * (M1**2*(math.sin(beta))**2 -1)/(M1**2*(g+math.cos(2*beta)) + 2) - math.tan(theta)
    beta, info, ier, msg = fsolve(equation, theta, full_output=True)
    if ier != 1:
        raise RuntimeError("No oblique shock solution")
    return beta[0]

def oblique_shock(theta: float, M1: float, g: float = DEFAULT_GAMMA):
    """
    Oblique Shock Relations
    
    :param theta: flow deflection angle [radians]
    :type theta: float
    :param M1: upstream Mach number
    :type M1: float
    :param g: gamma
    :type g: float
    :return: downstream Mach number, pressure ratio, density ratio, temperature ratio
    :rtype: tuple[float, float, float, float]
    """
    beta = oblique_shock_angle(theta, M1)
    M2_n, P2_P1, R2_R1, T2_T1 = normal_shock(M1*math.sin(beta))
    M2 = M2_n/(math.sin(beta-theta))
    return M2, P2_P1, R2_R1, T2_T1