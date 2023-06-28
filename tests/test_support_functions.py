from foam import support_functions as sf

def test_sign_pos():
    """ test sign positive"""
    result = sf.sign(1)
    expected = '+'
    assert result == expected

def test_sign_neg():
    """ test sign negative"""
    result = sf.sign(-5)
    expected = '-'
    assert result == expected

def test_get_param_from_filename():
    """ Test if the requested parameters get retrieved correctly from a filename."""
    file_path = 'Test-string_A7_B0_moreParams153_Comma3.56.hdf'
    parameters = ['A', 'B', 'Comma', 'absent']
    result = sf.get_param_from_filename(file_path, parameters)
    expected = {'A':'7', 'B':'0', 'Comma':'3.56'}
    assert result == expected

def test_split_line():
    line = 'string_to_be_split'
    result1, result2 = sf.split_line(line, 'to')
    expected1 = 'string_'
    expected2 = '_be_split'
    assert result1 == expected1
    assert result2 == expected2

def test_substring():
    line = 'string_to_be_split'
    result = sf.substring(line, 'st', '_')
    expected = 'ring'
    assert result == expected
