from compressible_flow.shocks import normal_shock

def test_normal_shock():
    M2, P2_P1, rho2_rho1, T2_T1 = normal_shock(2)
    assert M2 < 1
    assert P2_P1 > 1
    assert rho2_rho1 > 1
    assert T2_T1 > 1