from compressible_flow.airfoil import airfoil_cLcD

def test_airfoil_cL_cD():
    cl, cd = airfoil_cLcD(12, 2, 2)
    assert abs(cl) < 10
    assert abs(cd) < 10