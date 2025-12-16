import math
from .constants import DEFAULT_A, DEFAULT_GAMMA, DEFAULT_RHO
from .shocks import oblique_shock
from .expansions import expansion_fan

def freestream_pres(a: float = DEFAULT_A, rho: float = DEFAULT_RHO, g: float = DEFAULT_GAMMA) -> float:
    """
    Freestream pressure calculation
    
    :param a: speed of sound [m/s]
    :type a: float
    :param rho: density [kg/m^3]
    :type rho: float
    :param g: gamma
    :type g: float
    :return: freestream pressure [Pa]
    :rtype: float
    """
    return rho*a**2/g

def airfoil_pres(af_angle_deg: float, alpha_deg: float, M_inf: float,
                 a: float = DEFAULT_A, rho: float = DEFAULT_RHO, g: float = DEFAULT_GAMMA):
    """
    Diamond airfoil pressure coefficients on four panels (top: AB, BC; bottom: AD, DC)
    
    :param af_angle_deg: airfoil angle [degrees]
    :type af_angle_deg: float
    :param alpha_deg: angle of attack [degrees]
    :type alpha_deg: float
    :param M_inf: freestream Mach number
    :type M_inf: float
    :param a: speed of sound [m/s]
    :type a: float
    :param rho: density [kg/m^3]
    :type rho: float
    :param g: gamma
    :type g: float
    """
    theta = math.radians(af_angle_deg)
    alpha = math.radians(alpha_deg)
    v_inf = M_inf*a
    p_inf = freestream_pres(a, rho, g)
    
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
    else: #shock
        M_AB, P2_P1_AB, _, _ = oblique_shock(abs(delta_AB), M_inf)

    p_AB = p_inf * P2_P1_AB
    Cp_AB = (p_AB-p_inf)/(0.5*rho*v_inf**2)

    delta_BC = theta #always an expansion with positive AoA
    M_BC, P2_P1_BC, _, _ = expansion_fan(delta_BC, M_AB)
    p_BC = p_AB * P2_P1_BC
    Cp_BC = (p_BC-p_inf)/(0.5*rho*v_inf**2)
    return Cp_AB, Cp_BC, Cp_AD, Cp_DC

def airfoil_cLcD(af_angle_deg: float, alpha_deg: float, M_inf: float):
    """
    Compute CL and CD for a diamond-shaped airfoil
    
    :param af_angle_deg: airfoil angle [degrees]
    :type af_angle_deg: float
    :param alpha: angle of attack [degrees]
    :type alpha: float
    :param M_inf: freestream Mach number
    :type M_inf: float
    """
    cp_AB, cp_BC, cp_AD, cp_DC = airfoil_pres(af_angle_deg, alpha_deg, M_inf)

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