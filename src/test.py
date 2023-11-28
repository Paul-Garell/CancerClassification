import unittest
import seq


class TestAlign(unittest.TestCase):
    def test_basic_GATA3(self):
        GATA3_input1 = "CTGAAGAACATGAATTTCTCAAATCAAAGTAAAGCTGAGGAATCCTCTCCATGGCTCACTTCTACCTCCTCAATCTTAATGAGCTT"
        GATA3_input2 = (
            "GAGGTTGCAGTGAGCTGAGATCGCACCACTGCACTCCAGCCTGGGTGACAGAGCGAGACTCCATCT"
        )
        GATA3_input3 = "CCCTACTACGGAAACTCGGTCAGGGCCACGGTGCAGAGGTACCCTCCGACCCACCACGGTGAGTGCGCCCGGGGTGCCGGGGCTCCCGCCGGCCGCTTCAGCCGTCCCGGCTCGGGGAGGTCGGGAGG"

        self.assertEqual(seq.align(GATA3_input1), "GATA3")
        self.assertEqual(seq.align(GATA3_input2), "GATA3")
        self.assertEqual(seq.align(GATA3_input3), "GATA3")


# Running the tests
if __name__ == "__main__":
    unittest.main()
