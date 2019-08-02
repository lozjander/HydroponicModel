from nutrient_uptake_magnesium import mg_uptake

def test_mg_uptake():
    """First test for mg_uptake"""

    t = 999
    RSA = 999
    PAR = 999
    EC = 999
    a = 999
    b = 999
    c = 999
    d = 999
    DSR = 999
    uptake = mg_uptake(t, RSA, PAR, EC, a, b, c, d, DSR)

    assert uptake == 300


if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])


