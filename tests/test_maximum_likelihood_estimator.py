from foam import maximum_likelihood_estimator as mle
import numpy as np
import unittest

def test_merit_chi2():
    """ test chi2 merit function"""
    YObs  = np.asarray([1.0, 2.5, 3.7])
    ObsErr= np.asarray([0.1, 0.05, 0.3])
    YTheo = np.asarray([[1, 3, 4.6], [3, 2.1, 4.9], [1.0, 2.5, 4]])

    result = mle.merit_chi2(YObs, ObsErr, YTheo)
    expected = np.asarray([109, 480, 1])

    assert all ([abs((a-b)/a)< 1E-10 for a,b in zip(expected, result)])


def test_merit_mahalanobis():
    """ test MD merit function"""
    YObs  = np.asarray([1.0, 2.5, 3.7])
    ObsErr= np.asarray([0.1, 0.05, 0.3])
    YTheo = np.asarray([[1, 3, 4.6], [3, 2.1, 4.9], [1.0, 2.5, 4]])

    result = mle.merit_mahalanobis(YObs, ObsErr, YTheo, generate_output=False)
    expected = np.asarray([4.96270396, 4.997669, 0.91919192])

    assert all ([abs((a-b)/a)< 1E-8 for a,b in zip(expected, result)])        


class test_matrix(unittest.TestCase):
    def test_check_matrix_exit(self):
        matrix = np.asarray([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        with self.assertRaises(SystemExit) as cm:
            mle.check_matrix(matrix, generate_output=False)
        self.assertEqual(cm.exception.code, 1)

    def test_check_matrix(self):
        matrix = np.asarray([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
        mle.check_matrix(matrix, generate_output=False)
