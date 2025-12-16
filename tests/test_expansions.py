import math
from compressible_flow.expansions import expansion_fan
def test_expansion_fan():
    M2, P2_P1, rho2_rho1, T2_T1 = expansion_fan(math.radians(5), 2)
    assert M2 > 2
    assert P2_P1 < 1
    assert T2_T1 < 1