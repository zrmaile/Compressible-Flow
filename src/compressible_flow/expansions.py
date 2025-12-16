import math
from scipy.optimize import fsolve
from .constants import DEFAULT_GAMMA

def nu(M1: float, g: float = DEFAULT_GAMMA) -> float:
    """
    Prandtl-Meyer function nu(M) in radians
    
    :param M1: upstream Mach number
    :type M1: float
    :param g: gamma
    :type g: float
    :return: nu [radians]
    :rtype: float
    """
    return math.sqrt((g+1)/(g-1))*math.atan(math.sqrt((g-1)*(M1**2-1)/(g+1)))-math.atan(math.sqrt(M1**2-1))

def expansion_fan_M2(theta: float, M1:float, g: float = DEFAULT_GAMMA) -> float:
    """
    Solve for downstream Mach number, M2, given expansion angle theta and upstream Mach number
    
    :param theta: expansion angle
    :type theta: float
    :param M1: upstream Mach number
    :type M1: float
    :param g: gamma
    :type g: float
    :return: M2 - downstream Mach number
    :rtype: float
    """
    #Find nu(M2)
    nu_M2 = theta + nu(M1)
    #Find M2
    def equation(x):
        M = x[0]
        return math.sqrt((g+1)/(g-1))*math.atan(math.sqrt((g-1)*(M**2-1)/(g+1)))-math.atan(math.sqrt(M**2-1)) - nu_M2
    M2, info, ier, msg = fsolve(equation, 1.5, full_output=True)
    if ier != 1:
        raise RuntimeError("Error in isentropic flow equations")
    return M2[0]
def expansion_fan(theta: float, M1: float, g: float = DEFAULT_GAMMA):
    """
    Expansion fan relations
    
    :param theta: expansion angle
    :type theta: float
    :param M1: upstream Mach number
    :type M1: float
    :param g: gamma
    :type g: float
    :return: downstream Mach number, pressure ratio, density ratio, temperature ratio
    :rtype: tuple[float, float, float, float]
    """
    M2 = expansion_fan_M2(theta, M1)
    T2_T1 = (1 + (g-1)/2*M1**2)/(1 + (g-1)/2*M2**2)
    P2_P1 = (T2_T1)**(g/(g-1))
    R2_R1 = (P2_P1)**(1/g)
    return M2, P2_P1, R2_R1, T2_T1