from layer_photo import layer_photosyn


def test_layer_photo():
    """First test for layer_photo"""

    eps = 0.084
    P_gm = 2
    DIS_1 = 0.1127
    DIS_2 = 0.5
    DIS_3 = 0.8873
    WT_1 = 0.2778
    WT_2 = 0.4444
    WT_3 = 0.2778
    LAI = 2
    PAR = 200
    k_ext = 0.8
    photo = layer_photosyn(eps, P_gm, DIS_1, DIS_2, DIS_3, WT_1, WT_2, WT_3, LAI, PAR,
                           k_ext)

    assert photo == 3.690299845560409
    print(photo)


# #  Additional Code Added by Jens below
# T_average = 28.0
# dmp = dry_matter_production(W_p0, LAI0, leaf0, root0, stem0, fruit0, t, t_start,
#                             T_average, C_ppm, DIS_1, DIS_2, DIS_3, PAR, TEP,
#                             k_ext, WT_1, WT_2, WT_3, C_f, f)
#
# print("Dry matter production is: {}".format(dmp))

if __name__ == '__main__':
    import sys
    import pytest

    pytest.main([__file__] + sys.argv[1:])
