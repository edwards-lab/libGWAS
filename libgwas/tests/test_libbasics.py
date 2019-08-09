

import sys
import unittest

# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../../")
    sys.path.insert(0, "../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import libgwas

class TestBasics(unittest.TestCase):
    def testGenotypeData(self):
        gc = libgwas.GenotypeData()
        gc.append('0/0')
        gc.append('0/0')
        gc.append('0/1')
        gc.append('1/0')
        gc.append('0/0')
        gc.append('0/0')
        gc.append('1/1')

        self.assertEqual(7, len(gc.genotypes))
        self.assertEqual(2, gc.het_counts)
        self.assertEqual(10, gc.ref_counts)
        self.assertEqual(4, gc.alt_counts)
        self.assertEqual(0, gc.missing)

        self.assertAlmostEqual(4.0/14, gc.maf())
        self.assertAlmostEqual(4.0/14, gc.freq2())
        self.assertAlmostEqual(1-(4.0/14), gc.freq1())

    def testSysCall(self):
        cmd = f"wc -l {__file__}"
        wc = libgwas.sys_call(cmd)
        # I can see that there are more than 25 lines in the file
        self.assertTrue(int(wc.split()[0]) > 25)
        
        wc = libgwas.sys_call("wc -l non-existent-filename")
        self.assertEqual(None, wc)


if __name__ == "__main__":
    unittest.main()