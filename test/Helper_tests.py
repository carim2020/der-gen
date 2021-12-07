import unittest
import math
import numpy as np

from Helper import Vector3, cross_product, dot_product


class TestVector3(unittest.TestCase):
    def setUp(self):
        self.vector = Vector3()
        self.vector_x = Vector3(1, 0, 0)
        self.vector_y = Vector3(0, 1, 0)
        self.vector_z = Vector3(0, 0, 1)

    def takeDown(self):
        del self.vector
        del self.vector_x
        del self.vector_y
        del self.vector_z

    def test_str(self):
        self.assertEqual(str(self.vector), "0.0\t0.0\t0.0")
        self.assertEqual(str(self.vector_x), "1.0\t0.0\t0.0")
        self.assertEqual(str(self.vector_y), "0.0\t1.0\t0.0")
        self.assertEqual(str(self.vector_z), "0.0\t0.0\t1.0")

    def test_length_2(self):
        l1 = Vector3(1, 2, 0.5).length_2()
        self.assertEqual(l1, 5.25)
        self.assertEqual(self.vector.length_2(), 0)
        self.assertEqual(self.vector_x.length(), 1)
        self.assertEqual(self.vector_y.length(), 1)
        self.assertEqual(self.vector_z.length(), 1)
        self.assertNotEqual(l1, 5)

    def test_length(self):
        l1 = Vector3(1, 2, 0.5).length()
        self.assertAlmostEqual(l1, 2.2912878474779)
        self.assertEqual(self.vector.length_2(), 0)
        self.assertEqual(self.vector_x.length(), 1)
        self.assertEqual(self.vector_y.length(), 1)
        self.assertEqual(self.vector_z.length(), 1)
        self.assertNotEqual(l1, 5)

    def test_equal(self):
        self.assertTrue(Vector3(Vector3.MIN_DIST, 0, 0) == self.vector)
        self.assertTrue(Vector3(0, 1, 0) == self.vector_y)
        self.assertTrue(Vector3(0, 1 - Vector3.MIN_DIST, 0) == self.vector_y)
        self.assertTrue(Vector3(0, 1 + Vector3.MIN_DIST, 0) == self.vector_y)
        self.assertFalse(Vector3(0, 0, 0) == self.vector_y)

    def test_add(self):
        self.assertTrue(self.vector + self.vector_x == self.vector_x)
        self.assertTrue(self.vector + self.vector_y == self.vector_y)
        self.assertTrue(self.vector + self.vector_z == self.vector_z)
        self.assertFalse(self.vector + self.vector_y == self.vector_x)

    def test_sub(self):
        self.assertTrue(self.vector - self.vector_x == Vector3(-1, 0, 0))
        self.assertTrue(self.vector - self.vector_y == Vector3(0, -1, 0))
        self.assertTrue(self.vector - self.vector_z == Vector3(0, 0, -1))
        self.assertFalse(self.vector - self.vector_y == self.vector_x)

    def test_mul(self):
        self.assertTrue(self.vector_x * 5 == Vector3(5, 0, 0))
        self.assertTrue(self.vector_y * 4 == Vector3(0, 4, 0))
        self.assertTrue(self.vector_z * 2 == Vector3(0, 0, 2))
        self.assertTrue(Vector3(1, 2, 3) * 2 == Vector3(2, 4, 6))
        self.assertFalse(Vector3(0, 4, 3) * 3 == self.vector)

    def test_div(self):
        self.assertTrue(self.vector / 5 == self.vector)
        self.assertTrue(self.vector_y / 4 == Vector3(0, 0.25, 0))
        self.assertTrue(self.vector_z / 5 == Vector3(0, 0, 0.20))
        self.assertTrue(self.vector_x / 10 == Vector3(0.1, 0, 0))
        self.assertTrue(Vector3(1, 2, 3) / 2 == Vector3(0.5, 1, 1.5))
        self.assertFalse(Vector3(1, 2, 3) / 4 == Vector3(0.5, 0.5, 0.5))

    def test_iadd(self):
        v = Vector3(1, 1, 1)
        v += self.vector
        self.assertTrue(v == Vector3(1, 1, 1))
        self.assertFalse(v == self.vector)
        v += self.vector_x
        self.assertTrue(v == Vector3(2, 1, 1))

    def test_isub(self):
        v = Vector3(1, 1, 1)
        v -= self.vector
        self.assertTrue(v == Vector3(1, 1, 1))
        self.assertFalse(v == self.vector)
        v -= self.vector_x
        self.assertTrue(v == Vector3(0, 1, 1))

    def test_imul(self):
        v = Vector3(1, 1, 1)
        v *= 2
        self.assertTrue(v == Vector3(2, 2, 2))
        self.assertFalse(v == self.vector)

    def test_itruediv(self):
        v = Vector3(10, 10, 10)
        v /= 10
        self.assertTrue(v == Vector3(1, 1, 1))
        v /= 10
        self.assertFalse(v == Vector3(1, 1, 1))
        self.assertTrue(v == Vector3(0.1, 0.1, 0.1))

    def test_from_numpy(self):
        self.assertTrue(Vector3.from_numpy(np.array([0, 0, 0])) == self.vector)
        self.assertTrue(Vector3.from_numpy(np.array([1, 0, 0])) == self.vector_x)
        self.assertTrue(Vector3.from_numpy(np.array([0, 1, 0])) == self.vector_y)
        self.assertTrue(Vector3.from_numpy(np.array([0, 0, 1])) == self.vector_z)
        self.assertTrue(Vector3.from_numpy(np.array([1, 1, 1])) == Vector3(1, 1, 1))

    def test_to_numpy(self):
        self.assertListEqual(self.vector.to_numpy().tolist(), [0.0, 0.0, 0.0])
        self.assertListEqual(self.vector_x.to_numpy().tolist(), [1.0, 0.0, 0.0])
        self.assertListEqual(self.vector_y.to_numpy().tolist(), [0.0, 1.0, 0.0])
        self.assertListEqual(self.vector_z.to_numpy().tolist(), [0.0, 0.0, 1.0])
        self.assertListEqual(Vector3(1, 1, 1).to_numpy().tolist(), [1.0, 1.0, 1.0])

    def test_rotate_around_axis(self):
        v = Vector3(1.88, -2.85, 3.07)
        u = Vector3(3.04, 6.58, 4.26)
        v = v.rotate_around_axis(u, -math.pi/4)
        self.assertTrue(v == Vector3(-1.3902062920441292, -1.9027983900178758, 3.940619843692905))


class TestCrossAndDotProduct(unittest.TestCase):
    def setUp(self):
        self.a = Vector3(2, -1, 4)
        self.b = Vector3(-1, -2, -4)

    def test_cross_product(self):
        self.assertTrue(cross_product(self.a, self.b) == Vector3(12, 4, -5))

    def test_dot_product(self):
        self.assertTrue(dot_product(self.a, self.b) == -16)


if __name__ == '__main__':
    unittest.main()

