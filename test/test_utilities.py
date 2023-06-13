from unittest import TestCase, main, mock
import io
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


class Testget_version(TestCase):
    def _file_str(self, version):
        return "\n".join(("Docstring", 
                        "from . import predict",
                        "",
                        f"__version__ = {version}"))

    def test_version(self):
        for version in ["1.01.1dev1", "1.0.1", "0.0.0"]:

            with mock.patch("builtins.open", 
                            mock.mock_open(read_data=self._file_str(version))):

                self.assertEqual(version,
                        utils.get_version(file="irrelevant"))

    def test_incompatible_version(self):
        for version in ["01.1dev1", "1.0", "0"]:

            with mock.patch("builtins.open", 
                            mock.mock_open(read_data=self._file_str(version))):

                self.assertIsNone(utils.get_version(file="irrelevant"))
        




if __name__ == "__main__":
    main()
