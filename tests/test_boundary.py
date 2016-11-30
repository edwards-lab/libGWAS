#!/usr/bin/env python
import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../../")
    sys.path.insert(0, "../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")
from libgwas.boundary import BoundaryCheck
from libgwas.snp_boundary_check import SnpBoundaryCheck
import unittest


class TestBase(unittest.TestCase):
    def setUp(self):
        self.def_chrom = BoundaryCheck.chrom
        BoundaryCheck.chrom = -1     # Make sure it's set to an expected value


    def tearDown(self):
        BoundaryCheck.chrom = self.def_chrom

class TestBoundaryInitialization(TestBase):
    def testDefaultBoundaryInitialization(self):

        # By default, it will identify as invalid, since it didn't find any boundaries
        # This is just for simplifying command line parsing
        BoundaryCheck.chrom = -1
        b = BoundaryCheck()
        self.assertEqual(False, b.valid)

        # At this point, this should any valid chromosome/position combination
        self.assertTrue(b.TestBoundary(1, 100, ""))
        self.assertTrue(True, b.TestBoundary(10, 1000000, ""))
        self.assertTrue(True, b.TestBoundary(25, 10000, ""))

        # We should test that our short circuit functionality works
        self.assertTrue(b.NoExclusions())

    def testBoundaryInitSpaceAsSnp(self):

        # By default, it will identify as invalid, since it didn't find any boundaries
        # This is just for simplifying command line parsing
        BoundaryCheck.chrom = -1
        b = SnpBoundaryCheck(snps=[''])
        self.assertEqual(False, b.valid)

    def testBoundaryInitBP(self):
        BoundaryCheck.chrom = 1
        b = BoundaryCheck(bp=[10000, 500000])
        self.assertFalse(b.NoExclusions())
        self.assertTrue(b.valid)
        self.assertEqual(False, b.TestBoundary(1, 500, ""))
        self.assertEqual(True, b.TestBoundary(1, 10000, ""))
        self.assertEqual(True, b.TestBoundary(1, 500000, ""))
        self.assertEqual(True, b.TestBoundary(1, 250000, ""))
        self.assertEqual(False, b.TestBoundary(2, 250000, ""))
        self.assertEqual(False, b.TestBoundary(10, 10000, ""))

    def testBoundaryInitBPWithInclusions(self):
        BoundaryCheck.chrom = 1
        b = BoundaryCheck(bp=[10000, 500000])
        b.LoadSNPs(["rs12345", "rs23456"])
        self.assertFalse(b.NoExclusions())
        self.assertTrue(b.valid)
        self.assertEqual(False, b.TestBoundary(1, 500, ""))
        self.assertEqual(True, b.TestBoundary(1, 10000, ""))
        self.assertEqual(True, b.TestBoundary(1, 500000, ""))
        self.assertEqual(True, b.TestBoundary(1, 250000, ""))
        self.assertEqual(False, b.TestBoundary(2, 250000, ""))
        self.assertEqual(False, b.TestBoundary(10, 10000, ""))
        self.assertTrue(b.TestBoundary(1, 1000000, "rs12345"))
        self.assertTrue(b.TestBoundary(1, 1200000, "rs23456"))
        self.assertFalse(b.TestBoundary(1, 1200011, "rs345678"))

    def testBoundaryInitBPWithExclusions(self):
        BoundaryCheck.chrom = 1
        b = BoundaryCheck(bp=[10000, 500000])
        b.LoadExclusions(["rs12345", "rs234567", "rs345678"])
        self.assertFalse(b.NoExclusions())
        self.assertTrue(b.valid)
        self.assertFalse(b.TestBoundary(1, 500, ""))
        self.assertTrue(b.TestBoundary(1, 10000, ""))
        self.assertFalse(b.TestBoundary(1, 10010, "rs12345"))
        self.assertTrue(b.TestBoundary(1, 24000, "rs9876"))
        self.assertFalse(b.TestBoundary(1, 25000, "rs234567"))
        self.assertTrue(b.TestBoundary(1, 250000, ""))
        self.assertTrue(b.TestBoundary(1, 500000, ""))
        self.assertFalse(b.TestBoundary(2, 250000, ""))
        self.assertFalse(b.TestBoundary(10, 10000, ""))

    def testBoundaryInitKB(self):
        BoundaryCheck.chrom = 5
        b = BoundaryCheck(kb=[20, 50])
        self.assertFalse(b.NoExclusions())
        self.assertEqual(True, b.valid)
        self.assertEqual(False, b.TestBoundary(5, 15000, ""))
        self.assertEqual(True, b.TestBoundary(5, 20000, ""))
        self.assertEqual(True, b.TestBoundary(5, 30000, ""))
        self.assertEqual(True, b.TestBoundary(5, 50000, ""))
        self.assertEqual(False, b.TestBoundary(5, 50001, ""))
        self.assertEqual(False, b.TestBoundary(1, 25000, ""))
        self.assertEqual(False, b.TestBoundary(10, 20000, ""))

    def testBoundaryInitMB(self):
        BoundaryCheck.chrom = 10
        b = BoundaryCheck(mb=[1,3])
        self.assertTrue(b.valid)
        self.assertFalse(b.NoExclusions())
        self.assertTrue(b.TestBoundary(10, 1000000, ""))
        self.assertTrue(b.TestBoundary(10, 1200000, ""))
        self.assertTrue(b.TestBoundary(10, 3000000, ""))
        self.assertFalse(b.TestBoundary(10, 3000001, ""))
        self.assertFalse(b.TestBoundary(10, 999999, ""))
        self.assertFalse(b.TestBoundary(1, 1000500, ""))

    def testBoundaryExceedPos(self):
        BoundaryCheck.chrom = 10
        b = BoundaryCheck(mb=[1,3])
        self.assertTrue(b.valid)
        self.assertFalse(b.NoExclusions())
        self.assertFalse(b.TestBoundary(10, 100, ""))
        self.assertFalse(b.beyond_upper_bound)
        self.assertTrue(b.TestBoundary(10, 1000000, ""))
        self.assertFalse(b.beyond_upper_bound)
        self.assertTrue(b.TestBoundary(10, 1200000, ""))
        self.assertFalse(b.beyond_upper_bound)
        self.assertTrue(b.TestBoundary(10, 3000000, ""))
        self.assertFalse(b.beyond_upper_bound)
        self.assertFalse(b.TestBoundary(10, 3000001, ""))
        self.assertTrue(b.beyond_upper_bound)
        self.assertFalse(b.TestBoundary(10, 999999, ""))
        self.assertFalse(b.beyond_upper_bound)


class TestSnpBoundaryInitialization(TestBase):
    def testSnpRanges(self):
        BoundaryCheck.chrom = 22
        b = SnpBoundaryCheck(snps=["rs1-rs500","rs600-rs650","rs987654321"])
        self.assertFalse(b.NoExclusions())
        self.assertTrue(b.valid)
        self.assertFalse(b.TestBoundary(1, 10, "rs1"))
        self.assertFalse(b.TestBoundary(22, 20000, "rs750"))
        self.assertTrue(b.TestBoundary(22, 20001, "rs1"))
        self.assertTrue(b.TestBoundary(22, 22001, "rs50"))
        self.assertTrue(b.TestBoundary(22, 22005, "rs1"))
        # We don't really care which RS numbers we see, except for the boundaries
        self.assertTrue(b.TestBoundary(22, 22010, "rs1"))
        self.assertTrue(b.TestBoundary(22, 22011, "rs500"))
        self.assertFalse(b.TestBoundary(22, 23000, "rs499"))
        self.assertFalse(b.TestBoundary(22, 23001, "rs500"))
        self.assertTrue(b.TestBoundary(22, 23002, "rs600"))
        self.assertTrue(b.TestBoundary(22, 23003, "rs625"))
        self.assertTrue(b.TestBoundary(22, 23010, "rs20000"))
        self.assertTrue(b.TestBoundary(22, 24000, "rs650"))
        self.assertFalse(b.TestBoundary(22, 25000, "rs650"))
        self.assertFalse(b.TestBoundary(21, 2500000, "rs987654321"))
        self.assertTrue(b.TestBoundary(22, 2500001, "rs987654321"))

    def testSnpRangesWithExclusions(self):
        BoundaryCheck.chrom = 22
        b = SnpBoundaryCheck(snps=["rs1-rs500","rs600-rs650","rs987654321"])
        b.ignored_rs    = ["rs751","rs501"]
        self.assertFalse(b.NoExclusions())
        self.assertTrue(b.valid)
        self.assertFalse(b.TestBoundary(1, 10, "rs1"))
        self.assertFalse(b.TestBoundary(22, 20000, "rs750"))
        self.assertTrue(b.TestBoundary(22, 20001, "rs1"))
        self.assertFalse(b.TestBoundary(22, 20002, "rs751"))
        self.assertTrue(b.TestBoundary(22, 22003, "rs50"))
        self.assertTrue(b.TestBoundary(22, 22005, "rs1"))
        # We don't really care which RS numbers we see, except for the boundaries
        self.assertTrue(b.TestBoundary(22, 22010, "rs1"))
        self.assertFalse(b.TestBoundary(22, 22011, "rs501"))
        self.assertTrue(b.TestBoundary(22, 22012, "rs500"))
        self.assertFalse(b.TestBoundary(22, 23000, "rs499"))
        self.assertFalse(b.TestBoundary(22, 23002, "rs500"))
        self.assertTrue(b.TestBoundary(22, 23003, "rs600"))
        self.assertTrue(b.TestBoundary(22, 23004, "rs625"))
        self.assertTrue(b.TestBoundary(22, 23010, "rs20000"))
        self.assertTrue(b.TestBoundary(22, 24000, "rs650"))
        self.assertFalse(b.TestBoundary(22, 25000, "rs650"))
        self.assertFalse(b.TestBoundary(21, 2500000, "rs987654321"))
        self.assertTrue(b.TestBoundary(22, 2500000, "rs987654321"))

class TestSnpBoundaryInitializationNoChr(TestBase):
    """We are labelling the missing chromosomes in MACH as NA (as are positions).
    """
    def testSnpRanges(self):
        BoundaryCheck.chrom = "NA"
        b = SnpBoundaryCheck(snps=["rs1-rs500","rs600-rs650","rs987654321"])
        self.assertFalse(b.NoExclusions())
        self.assertTrue(b.valid)
        self.assertFalse(b.TestBoundary(1, 10, "rs1"))
        self.assertFalse(b.TestBoundary("NA", "NA", "rs750"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs1"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs50"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs1"))
        # We don't really care which RS numbers we see, except for the boundaries
        self.assertTrue(b.TestBoundary("NA", "NA", "rs1"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs500"))
        self.assertFalse(b.TestBoundary("NA", "NA", "rs499"))
        self.assertFalse(b.TestBoundary(22, 23001, "rs500"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs600"))
        self.assertTrue(b.TestBoundary("NA", '23003', "rs625"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs20000"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs650"))
        self.assertFalse(b.TestBoundary("NA", "NA", "rs650"))
        self.assertFalse(b.TestBoundary(21, "NA", "rs987654321"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs987654321"))

    def testSnpRangesWithExclusions(self):
        BoundaryCheck.chrom = "NA"
        b = SnpBoundaryCheck(snps=["rs1-rs500","1:600-1:650","1:987654321"])
        b.ignored_rs    = ["1:751","1:501"]
        self.assertFalse(b.NoExclusions())
        self.assertTrue(b.valid)
        self.assertFalse(b.TestBoundary(1, "NA", "rs1"))
        self.assertFalse(b.TestBoundary("NA", "NA", "rs750"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs1"))
        self.assertFalse(b.TestBoundary("NA", "NA", "1:751"))
        self.assertTrue(b.TestBoundary("NA", "NA", "1:50"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs1"))
        # We don't really care which RS numbers we see, except for the boundaries
        self.assertTrue(b.TestBoundary("NA", "NA", "rs1"))
        self.assertFalse(b.TestBoundary("NA", "NA", "1:501"))
        self.assertTrue(b.TestBoundary("NA", "NA", "rs500"))
        self.assertFalse(b.TestBoundary("NA", "NA", "rs499"))
        self.assertFalse(b.TestBoundary("NA", "NA", "rs500"))
        self.assertTrue(b.TestBoundary("NA", "NA", "1:600"))
        self.assertTrue(b.TestBoundary("NA", "NA", "1:625"))
        self.assertTrue(b.TestBoundary("NA", "NA", "1:20000"))
        self.assertTrue(b.TestBoundary("NA", "NA", "1:650"))
        self.assertFalse(b.TestBoundary("NA", "NA", "1:650"))
        self.assertFalse(b.TestBoundary(21, "NA", "1:987654321"))
        self.assertTrue(b.TestBoundary("NA", "NA", "1:987654321"))




if __name__ == "__main__":
    unittest.main()
