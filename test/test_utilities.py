from unittest import TestCase, main
import numpy as np
from afcn import utils


class TestIsNumericNparray(TestCase):
    def test_true_numeric(self):
        self.assertTrue(utils.is_numeric_nparray(np.array([1,2])))
        self.assertTrue(utils.is_numeric_nparray(np.array([1,False])))

    def test_False_numeric(self):
        self.assertFalse(utils.is_numeric_nparray(np.array([False])))
        self.assertFalse(utils.is_numeric_nparray(np.array(["1",2])))
        self.assertFalse(utils.is_numeric_nparray(np.array(["1",2.1])))
        self.assertFalse(utils.is_numeric_nparray(np.array(["1",np.pi])))
        self.assertFalse(utils.is_numeric_nparray(np.array([None,2])))
        self.assertFalse(utils.is_numeric_nparray(np.array([None,False])))
        self.assertFalse(utils.is_numeric_nparray(np.array([True,False])))



class TestIsBiallelic(TestCase):
    def test_true_cases(self):
        self.assertTrue(utils.is_biallelic(np.array([1,1,1,1])))
        self.assertTrue(utils.is_biallelic(np.array([1.,1.,1.,1.])))
        self.assertTrue(utils.is_biallelic(np.array([0,0,0,0])))
        self.assertTrue(utils.is_biallelic(np.array([1,0,0,0])))
        self.assertTrue(utils.is_biallelic(np.array([1,0,0,1])))
        self.assertTrue(utils.is_biallelic([1,0,0,1]))
        self.assertTrue(utils.is_biallelic(np.array([1])))
        self.assertTrue(utils.is_biallelic(np.array([0])))
        self.assertTrue(utils.is_biallelic(0))
        self.assertTrue(utils.is_biallelic(1))
        self.assertTrue(utils.is_biallelic(np.array([0,1,True])))

    def test_false_cases(self):
        self.assertFalse(utils.is_biallelic([0, 1, 2]))
        self.assertFalse(utils.is_biallelic(np.array([0,1,2])))
        self.assertFalse(utils.is_biallelic(np.array([0,1,1.1])))
        self.assertFalse(utils.is_biallelic(np.array([0,1,1, np.nan])))
        self.assertFalse(utils.is_biallelic(np.array([0,1,1, np.nan])))
        self.assertFalse(utils.is_biallelic(2))
        self.assertFalse(utils.is_biallelic(np.nan))
        self.assertFalse(utils.is_biallelic(None))
        self.assertFalse(utils.is_biallelic(True))
        self.assertFalse(utils.is_biallelic(False))
        self.assertFalse(utils.is_biallelic(np.array(["1", "0"])))

    def test_exception(self):
        with self.assertRaises(TypeError):
            utils.is_biallelic(np.array([0,None,2]))


if __name__ == "__main__":
    main()
