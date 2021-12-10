import unittest


from Translator import ob2dergen
from openbabel import openbabel as ob


class TestSubTemplate(unittest.TestCase):
    def setUp(self) -> None:
        conv = ob.OBConversion()
        conv.SetInFormat("xyz")
        obm = ob.OBMol()
        conv.ReadFile(obm, "../res/2-propylylcyclopentenone.xyz")
        self.mol = ob2dergen(obm)
        self.unit_v = (self.mol.get_atom(6).coord - self.mol.get_atom(1).coord).unit()

    def test_init(self):
        raise NotImplementedError
        # red = Nitrile(self.mol)
        # self.assertIsInstance(red, Nitrile)

    def test_generate_new_compound(self):
        raise NotImplementedError
        # red = Nitrile(self.mol)
        # red.generate_new_compound(1, 0, 6)
        #
        # self.assertEqual(len(red.atoms), 9)
        # self.assertEqual(red.get_atom(8), Atom(symbol=ATOMIC_SYMBOLS[6], atomic_number=6,
        #                                        coord=Vector3(-1.9884643835538665, -0.8483596826630728,
        #                                                      0.9720840138877235)))
        # self.assertEqual(red.get_atom(9), Atom(symbol=ATOMIC_SYMBOLS[7], atomic_number=7,
        #                                        coord=Vector3(-2.873585852121078, -0.567382946071825,
        #                                                      1.633304143175962)))
        #
        # self.assertEqual(red.bonds["8_9"], BondType.TRIPLE)
        # self.assertEqual(red.bonds["1_8"], BondType.SINGLE)


if __name__ == "__main__":
    unittest.main()
