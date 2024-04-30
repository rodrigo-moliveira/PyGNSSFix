import unittest


def a(x):
    # Implementation of function a
    return x * 2


def b(x):
    # Implementation of function b
    return x + x


class TestFunctions(unittest.TestCase):

    def test_output(self):
        # Define input
        input_value = 5

        # Define threshold
        threshold = 10

        # Get output of functions a and b
        output_a = a(input_value)
        output_b = b(input_value)

        # Assert that outputs are equal within the threshold
        self.assertAlmostEqual(output_a, output_b, delta=threshold)


if __name__ == '__main__':
    unittest.main()
