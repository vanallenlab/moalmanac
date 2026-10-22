import unittest

from illustrator import PreclinicalEfficacy


class UnitTestPreclinicalEfficacy(unittest.TestCase):
    def test_create_y_range(self):
        observed = PreclinicalEfficacy.create_y_range([0.05, 3.2, 150.0])
        self.assertEqual(observed, [0.01, 0.1, 1, 10, 100, 1000])

    def test_create_y_range_ignores_non_positive(self):
        observed = PreclinicalEfficacy.create_y_range([0.0, 0.05, 3.2])
        self.assertEqual(observed, [0.01, 0.1, 1, 10])

    def test_create_y_range_no_positive_values(self):
        self.assertEqual(PreclinicalEfficacy.create_y_range([0.0, 0.0]), [])
        self.assertEqual(PreclinicalEfficacy.create_y_range([]), [])


if __name__ == "__main__":
    unittest.main()
