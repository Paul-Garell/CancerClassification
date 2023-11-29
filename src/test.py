import unittest
import seq

# Run Using python src/test.py


class TestAlign(unittest.TestCase):
    def test_identity_GATA3(self):
        GATA3_input1 = "CTGAAGAACATGAATTTCTCAAATCAAAGTAAAGCTGAGGAATCCTCTCCATGGCTCACTTCTACCTCCTCAATCTTAATGAGCTT"
        GATA3_input2 = (
            "GAGGTTGCAGTGAGCTGAGATCGCACCACTGCACTCCAGCCTGGGTGACAGAGCGAGACTCCATCT"
        )
        GATA3_input3 = "CCCTACTACGGAAACTCGGTCAGGGCCACGGTGCAGAGGTACCCTCCGACCCACCACGGTGAGTGCGCCCGGGGTGCCGGGGCTCCCGCCGGCCGCTTCAGCCGTCCCGGCTCGGGGAGGTCGGGAGG"

        self.assertEqual(seq.align(GATA3_input1), ["GATA3", 8113468, 8113554, []])
        self.assertEqual(seq.align(GATA3_input2), ["GATA3", 8110649, 8110715, []])
        self.assertEqual(seq.align(GATA3_input3), ["GATA3", 8097802, 8097930, []])

    def test_change_GATA3(self):
        GATA3_input1 = "CAGAAGAACATGAATTTCTCAAATCAAAGTAAAGCTGAGGAATCCTCTCCATGGCTCACTTCTACCTCCTCAATCTTAATGAGCTT"

        self.assertEqual(
            seq.align(GATA3_input1), ["GATA3", 8113468, 8113554, [("T", "A")]]
        )


# Running the tests
if __name__ == "__main__":
    unittest.main()
