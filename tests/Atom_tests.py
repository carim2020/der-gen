import unittest

from Helper import Vector3
from Atom import Atom


class TestAtom(unittest.TestCase):
    def setUp(self):
        self.atom = Atom()
        self.atom.coord = Vector3(1, 1, 1)
        self.atom.symbol = "C"
        self.atom.atomic_num = 6
        self.atom.atomic_mass = 12.0
        self.atom.id = 2

    def takeDown(self):
        del self.atom

    def test_init(self) -> None:
        self.assertIsInstance(self.atom, Atom)

    def test_str(self) -> None:
        self.assertEqual(str(self.atom), f"C\t{Vector3(1, 1, 1)}")

    def test_coord(self) -> None:
        self.assertEqual(self.atom.coord, Vector3(1, 1, 1))
        self.atom.coord = Vector3(2, 2, 2)
        self.assertNotEqual(self.atom.coord, Vector3(1, 1, 1))
        self.assertEqual(self.atom.coord, Vector3(2, 2, 2))
        self.atom.coord = Vector3(1, 1, 1)

    def test_symbol(self) -> None:
        self.assertEqual(self.atom.symbol, "C")
        self.atom.symbol = "N"
        self.assertNotEqual(self.atom.symbol, "C")
        self.assertEqual(self.atom.symbol, "N")
        self.atom.symbol = "C"

    def test_atomic_num(self) -> None:
        self.assertIsInstance(self.atom.atomic_num, int)
        self.assertEqual(self.atom.atomic_num, 6)
        self.atom.atomic_num = 7
        self.assertNotEqual(self.atom.atomic_num, 6)
        self.assertEqual(self.atom.atomic_num, 7)

    def test_atomic_mass(self) -> None:
        self.assertIsInstance(self.atom.atomic_mass, float)
        self.assertEqual(self.atom.atomic_mass, 12.0)
        self.atom.atomic_mass = 14.0
        self.assertNotEqual(self.atom.atomic_mass, 12.0)
        self.assertEqual(self.atom.atomic_mass, 14.0)
        self.atom.atomic_mass = 12.0

    def test_id(self) -> None:
        self.assertIsInstance(self.atom.id, int)
        self.assertEqual(self.atom.id, 2)
        self.atom.id = 3
        self.assertNotEqual(self.atom.id, 2)
        self.assertEqual(self.atom.id, 3)
        self.atom.id = 2


if __name__ == '__main__':
    unittest.main()
