import numpy as np
from foam import build_optimised_pattern as bop

def test_generate_spacing_series():
    """Test period spacing series calculated by generate_spacing_series"""
    periods = [0.2, 0.3, 0.38, 0.42, 0.45]
    expected = ([0.1*86400, 0.08*86400, 0.04*86400, 0.03*86400], None)
    result = bop.generate_spacing_series(periods)

    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected[0], result[0])])
    assert result[1] == expected[1]

def test_generate_spacing_series_with_errors():
    """Test period spacing series and its errors calculated by test_generate_spacing_series"""
    periods = [1.1, 1.3, 1.39, 1.47, 1.52 ]
    errors = [1E-1, 2E-1, 4E-1, 3E-1, 2E-1]
    expected_spacing = [0.2*86400, 0.09*86400, 0.08*86400, 0.05*86400]
    expected_errors  = [np.sqrt(0.05)*86400,  np.sqrt(0.2)*86400, np.sqrt(0.25)*86400, np.sqrt(0.13)*86400]
    result = bop.generate_spacing_series(periods, errors)

    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected_spacing, result[0])])
    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected_errors, result[1])])

def test_puls_series_from_given_puls_at_end():
    """Test building theoretical pulsation pattern from end of pattern."""
    theory_freq = np.asarray([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
    obs_freq  = np.asarray([0.68, 0.76, 0.88, 0.95, 1.02, 1.06])
    obs_to_build_from = 1.06
    expected = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
    result = bop.puls_series_from_given_puls(theory_freq, obs_freq, obs_to_build_from)

    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected, result)])

def test_puls_series_from_given_puls_midway():
    """Test building theoretical pulsation pattern from midway pattern."""
    theory_freq = np.asarray([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
    obs_freq = np.asarray([0.68, 0.76, 0.88, 0.95, 1.02, 1.06])
    obs_to_build_from = 0.88
    expected = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
    result = bop.puls_series_from_given_puls(theory_freq, obs_freq, obs_to_build_from)
    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected, result)])

def test_chisq_longest_sequence():
    """Test building theoretical pulsation pattern according to the longest sequence method."""
    theory_freq = np.asarray([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
    orders    = np.asarray([-1, -2, -3, -4, -5, -6, -7, -8])
    obs_freq  = np.asarray([0.68, 0.76, 0.88, 0.95, 1.02, 1.06])
    obs_errs  = np.asarray([ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
    expected_freq = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
    expected_orders = [-2, -3, -4, -5, -6, -7]
    result = bop.chisq_longest_sequence(theory_freq, orders, obs_freq, obs_errs)
    assert abs(result[0]-6.199999999999989) < 1E-15
    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected_freq, result[1])])
    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected_orders, result[2])])

def test_chisq_longest_sequence2():
    """Test building theoretical pulsation pattern according to the longest sequence method."""
    theory_freq = np.asarray([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
    orders    = np.asarray([-1, -2, -3, -4, -5, -6, -7, -8])
    obs_freq  = np.asarray([0.68, 0.76, 0.84, 0.94, 1.02, 1.06])
    obs_errs  = np.asarray([ 0.05, 0.03, 0.07, 0.08, 0.06, 0.1])
    expected_freq = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
    expected_orders = [-1, -2, -3, -4, -5, -6]
    result = bop.chisq_longest_sequence(theory_freq, orders, obs_freq, obs_errs)
    assert abs(result[0]-0.09826369168357019) < 1E-17
    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected_freq, result[1])])
    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected_orders, result[2])])    
