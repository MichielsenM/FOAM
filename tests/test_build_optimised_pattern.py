from foam import build_optimised_pattern as bop

def test_generate_spacing_series():
    """ test period spacing series calculated by generate_spacing_series"""
    periods = [0.2, 0.3, 0.38, 0.42, 0.45]
    expected = ([0.1*86400, 0.08*86400, 0.04*86400, 0.03*86400], None)
    result = bop.generate_spacing_series(periods)

    assert all ([(a-b)/a< 1E-10 for a,b in zip(expected[0], result[0])])
    assert result[1] == expected[1]

def test_generate_spacing_series_with_errors():
    """ test period spacing series and its errors calulated by test_generate_spacing_series """
    periods = [1.1, 1.3, 1.39, 1.47, 1.52 ]
    errors = [1E-6, 2E-7, 6E-6, 4E-6, 2E-6]
    expected = ([0.2*86400, 0.09*86400, 0.08*86400, 0.05*86400], [8.9856E-8*86400, 3.113856E-6*86400, 4.4928E-6*86400, 1.728E-8*86400])
    result = bop.generate_spacing_series(periods, errors)

    assert all ([(a-b)/a< 1E-10 for a,b in zip(expected[0], result[0])])
    assert all ([(a-b)/a< 1E-10 for a,b in zip(expected[1], result[1])])
