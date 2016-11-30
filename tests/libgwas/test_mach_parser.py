#!/usr/bin/env python

import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../../")
    sys.path.insert(0, "../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")


import unittest
import numpy
import os

from libgwas.data_parser import DataParser
from libgwas.pheno_covar import PhenoCovar
from libgwas import mach_parser
from libgwas.boundary import BoundaryCheck


import gzip

base_freq = [0.95, 0.75, 0.7,0.8, 0.65, 0.7, 0.85, 0.7, 0.7, 0.3]

numpy.random.seed(1557)

class TestBase(unittest.TestCase):
    def setUp(self):
        self.allele_1 = list("AAACCCGGGTCGTGTATACC")
        self.allele_2 = list("CGTGTATACCAAACCCGGGT")

        self.WriteTestFiles()

        self.phenotypes = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5, 0.9, 1.0]
        self.sex = [1,2,1,1,2,2,1,1,2,2,2,1]
        self.chrpos_encoding = mach_parser.Parser.chrpos_encoding
        self.dosage_ext = mach_parser.Parser.dosage_ext
        self.info_ext = mach_parser.Parser.info_ext
        self.chrom = BoundaryCheck.chrom
        self.boundary = DataParser.boundary
        DataParser.boundary = BoundaryCheck()
        self.min_maf = DataParser.min_maf
        self.max_maf = DataParser.max_maf
        self.snp_miss_tol = DataParser.snp_miss_tol
        self.ind_miss_tol = DataParser.ind_miss_tol
        DataParser.ind_exclusions = []
        self.sex_as_covar = PhenoCovar.sex_as_covariate
        self.has_sex = DataParser.has_sex
        self.has_pheno = DataParser.has_pheno
        self.has_parents = DataParser.has_parents
        self.has_fid = DataParser.has_fid
        self.has_liability = DataParser.has_liability
        self.encoding = mach_parser.encoding
        self.compression = DataParser.compressed_pedigree
        DataParser.compressed_pedigree = True
        #self.chunk_stride = mach_parser.chunk_stride


    def tearDown(self):
        os.remove(self.gen_file)
        os.remove(self.gen_file2)
        os.remove(self.uncmp_1)
        os.remove(self.uncmp_2)
        os.remove(self.info_file1)
        os.remove(self.info_file2)
        os.remove(self.info_ucmp1)
        os.remove(self.info_ucmp2)

        mach_parser.Parser.dosage_ext = self.dosage_ext
        mach_parser.Parser.info_ext = self.info_ext
        mach_parser.Parser.chrpos_encoding = self.chrpos_encoding

        BoundaryCheck.chrom  = self.chrom
        DataParser.boundary  = self.boundary
        DataParser.min_maf   = self.min_maf
        DataParser.max_maf   = self.max_maf
        DataParser.snp_miss_tol  = self.snp_miss_tol
        DataParser.ind_miss_tol  = self.ind_miss_tol
        DataParser.ind_exclusions = []
        DataParser.has_sex    = self.has_sex
        DataParser.has_pheno  = self.has_pheno
        DataParser.has_fid    = self.has_fid
        DataParser.has_liability = self.has_liability
        DataParser.has_parents = self.has_parents
        PhenoCovar.sex_as_covariate = self.sex_as_covar
        DataParser.compressed_pedigree = self.compression
        #mach_parser.chunk_stride = self.chunk_stride


    def WriteTestFiles(self, prefix = "__test_imputed"):

        self.ind_ids = ["ID0001:FAM001",
                        "ID0002:FAM002",
                        "ID0003:FAM003",
                        "ID0004:FAM004",
                        "ID0005:FAM005",
                        "ID0006:FAM006",
                        "ID0007:FAM007",
                        "ID0008:FAM008",
                        "ID0009:FAM009",
                        "ID0010:FAM010",
                        "ID0011:FAM011",
                        "ID0012:FAM012"]

        self.gen_file = "%s.dose.gz" % (prefix)
        self.gen_file2 = "%s-2.dose.gz" % (prefix)
        self.info_file1 = "%s.info.gz" % (prefix)
        self.info_file2 = "%s-2.info.gz" % (prefix)
        self.uncmp_1 = "%s.dose" % (prefix)
        self.uncmp_2 = "%s-2.dose" % (prefix)
        self.info_ucmp1 = "%s.info" % (prefix)
        self.info_ucmp2 = "%s-2.info" % (prefix)
        gen_file = gzip.open(self.gen_file, 'wb')
        uncmp_file = open(self.uncmp_1, 'w')
        idx = 0
        self.dosage_encoding = numpy.zeros((20, 12))
        self.positions = []
        self.mafs = numpy.zeros(len(base_freq) * 2)

        info_file = open(self.info_file1, 'w')
        print >> info_file, "snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0"

        self.chroms = [ int(x) for x in ['1'] * 7 + ['2'] * 7 + ['3'] * 6]
        self.positions = [1012, 1020, 1026, 1032, 1100, 1137, 1149] * 2 + [1012, 1020, 1026, 1032, 1100, 1137]
        self.alleles = [list(numpy.random.choice(['A','C','G','T'], 2, replace=False)) for x in range(0, 20)]
        idx = 0

        mafs = numpy.zeros((10))
        dosages = numpy.zeros((12, 10))
        for ind in self.ind_ids:
            f = numpy.random.normal(base_freq, scale=0.1)
            f[f>=1.0] = 0.99
            maf = 1.0 - f
            AA = f * f
            Aa = 2 * f * maf
            aa = maf * maf
            dosages[idx] = Aa + 2*AA
            mafs += dosages[idx] / 2
            print >> gen_file, "\t".join([
                ind,
                "DOSE"] +
                ["%.3f" % x for x in dosages[idx]]
            )
            print >> uncmp_file, "\t".join([
                ind,
                "DOSE"] +
                ["%.3f" % x for x in dosages[idx]]
            )
            idx += 1
        self.mafs[0:10] = mafs/10
        self.dosage_encoding[0:10,:] = numpy.transpose(dosages)
        gen_file.close()
        uncmp_file.close()
        info_file = gzip.open(self.info_file1, 'wb')
        info_ufile = open(self.info_ucmp1, 'w')
        print >> info_file, "SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tdose2"
        print >> info_ufile, "SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tdose2"
        for idx in range(0, 10):
            print >> info_file, "\t".join([
                "%s:%d" % (self.chroms[idx],self.positions[idx]),
                self.allele_1[idx],
                self.allele_2[idx],
                str(1.0-self.mafs[idx]),
                str(self.mafs[idx]),
                '0.99912',
                '0.8',
                "\t".join(['-'] * 6)
            ])
            print >> info_ufile, "\t".join([
                "%s:%d" % (self.chroms[idx],self.positions[idx]),
                self.allele_1[idx],
                self.allele_2[idx],
                str(1.0-self.mafs[idx]),
                str(self.mafs[idx]),
                '0.99912',
                '0.8',
                "\t".join(['-'] * 6)
            ])
        info_file.close()
        info_ufile.close()


        gen_file = gzip.open(self.gen_file2, 'wb')
        uncmp_file = open(self.uncmp_2, 'w')

        idx = 0
        mafs = numpy.zeros((10))
        dosages = numpy.zeros((12, 10))
        for ind in self.ind_ids:
            f = numpy.random.normal(base_freq, scale=0.1)
            f[f>=1.0] = 0.99
            maf = 1.0 - f
            mafs += maf
            AA = f * f
            Aa = 2 * f * maf
            aa = maf * maf
            dosages[idx] = Aa + 2*aa
            print >> gen_file, "\t".join([
                ind,
                "DOSE"] +
                ["%.3f" % x for x in dosages[idx]]
            )
            print >> uncmp_file, "\t".join([
                ind,
                "DOSE"] +
                ["%.3f" % x for x in dosages[idx]]
            )
            idx += 1
        self.mafs[10:] = mafs/10
        self.dosage_encoding[10:,:] = numpy.transpose(dosages)

        gen_file.close()

        info_file = gzip.open(self.info_file2, 'wb')
        info_cfile = open(self.info_ucmp2, 'w')
        print >> info_file, "SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tdose2"
        print >> info_cfile, "SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tdose2"
        for idx in range(10, 20):
            print >> info_file, "\t".join([
                "%s:%d" % (self.chroms[idx],self.positions[idx]),
                self.allele_1[idx],
                self.allele_2[idx],
                str(1.0-self.mafs[idx]),
                str(self.mafs[idx]),
                '0.99912',
                '0.8',
                "\t".join(['-'] * 6)
            ])
            print >> info_cfile, "\t".join([
                "%s:%d" % (self.chroms[idx],self.positions[idx]),
                self.allele_1[idx],
                self.allele_2[idx],
                str(1.0-self.mafs[idx]),
                str(self.mafs[idx]),
                '0.99912',
                '0.8',
                "\t".join(['-'] * 6)
            ])
        info_cfile.close()
        info_file.close()

        self.mach_parser = mach_parser.Parser([self.gen_file])

class TestImputedBasics(TestBase):
    def testChromosomes(self):
        mach_parser.Parser.chrpos_encoding = True

        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(self.chroms[idx], snp.chr)

            idx += 1
        self.assertEqual(20, idx)
    def testChromosomesNoChrPos(self):
        mach_parser.Parser.chrpos_encoding = False

        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual("NA", snp.pos)
            self.assertEqual("NA", snp.chr)
            self.assertEqual("%d:%d" % (self.chroms[idx], self.positions[idx]), snp.rsid)
            idx += 1
        self.assertEqual(20, idx)
    def testValuesUncompressed(self):
        mach_parser.Parser.chrpos_encoding = True
        DataParser.compressed_pedigree = False
        mach_parser.Parser.gen_ext = "gen"
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.uncmp_1, self.uncmp_2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(self.chroms[idx], snp.chr)

            for i in range(0, len(self.dosage_encoding[idx])):
                self.assertAlmostEqual(self.dosage_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)


    def testValues(self):
        mach_parser.Parser.chrpos_encoding = True
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(self.chroms[idx], snp.chr)

            for i in range(0, len(self.dosage_encoding[idx])):
                self.assertAlmostEqual(self.dosage_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)


    def testInfoFileUse(self):
        mach_parser.Parser.chrpos_encoding = True

        # We'll give it an invalid gen_ext so that we can be certain that it's using the files provided
        mach_parser.Parser.gen_ext='asdf'
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2], info_files=[self.info_file1, self.info_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.dosage_encoding[idx])):
                self.assertAlmostEqual(self.dosage_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testInfoFileUseNoChrPos(self):
        # We'll give it an invalid gen_ext so that we can be certain that it's using the files provided
        mach_parser.Parser.gen_ext='asdf'
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2], info_files=[self.info_file1, self.info_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual("NA", snp.pos)
            self.assertEqual("NA", snp.chr)
            self.assertEqual("%s:%s" % (self.chroms[idx], self.positions[idx]), snp.rsid)
            for i in range(0, len(self.dosage_encoding[idx])):
                self.assertAlmostEqual(self.dosage_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testMAF(self):
        mach_parser.Parser.chrpos_encoding = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            maf = numpy.mean(snp.genotype_data/2)
            self.assertAlmostEqual(maf, snp.maf, places=3)
            idx += 1
        self.assertEqual(10, idx)

    def testAlelles(self):
        mach_parser.Parser.chrpos_encoding = True

        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(snp.major_allele, self.allele_2[idx])
            self.assertEqual(snp.minor_allele, self.allele_1[idx])
            idx += 1
        self.assertEqual(10, idx)

    def testChunkStride(self):
        mach_parser.Parser.chrpos_encoding = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file])
        parser.chunk_stride = 7


        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(snp.major_allele, self.allele_2[idx])
            self.assertEqual(snp.minor_allele, self.allele_1[idx])
            idx += 1
        self.assertEqual(10, idx)
    def testLongerList(self):
        mach_parser.Parser.chrpos_encoding = True
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2]*3)
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        values = numpy.vstack((self.dosage_encoding, self.dosage_encoding, self.dosage_encoding))
        positions = self.positions * 3

        for snp in parser:
            self.assertEqual(positions[idx], snp.pos)
            for i in range(0, len(values[idx])):
                self.assertAlmostEqual(values[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(60, idx)

    def testFilteredInd(self):
        mach_parser.Parser.chrpos_encoding = True
        DataParser.ind_exclusions = self.ind_ids[0:2]
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(10, snp.genotype_data.shape[0])
            for i in range(2, len(self.dosage_encoding[idx])):
                self.assertAlmostEqual(self.dosage_encoding[idx][i], snp.genotype_data[i-2], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testBoundaried(self):
        mach_parser.Parser.chrpos_encoding = True

        BoundaryCheck.chrom = 1
        DataParser.boundary = BoundaryCheck(bp=[0, 1137])
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.dosage_encoding[idx])):
                self.assertAlmostEqual(self.dosage_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(6, idx)

    def testBoundedUpper(self):
        mach_parser.Parser.chrpos_encoding = True

        BoundaryCheck.chrom = 3
        DataParser.boundary = BoundaryCheck(bp=[1026, 2000])
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 16

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.dosage_encoding[idx])):
                self.assertAlmostEqual(self.dosage_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testBoundedMiddle(self):
        mach_parser.Parser.chrpos_encoding = True
        BoundaryCheck.chrom = 2
        DataParser.boundary = BoundaryCheck(bp=[1020, 1137])
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file, self.gen_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 8
        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.dosage_encoding[idx])):
                self.assertAlmostEqual(self.dosage_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(13, idx)


    def testFilterMAF(self):
        mach_parser.Parser.chrpos_encoding = True
        DataParser.min_maf = 0.45
        pc = PhenoCovar()
        parser = mach_parser.Parser([self.gen_file])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0
        for snp in parser:
            while numpy.mean(self.mafs[idx]) < DataParser.min_maf:
                idx += 1
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(snp.major_allele, self.allele_2[idx])
            self.assertEqual(snp.minor_allele, self.allele_1[idx])
            idx += 1




if __name__ == "__main__":
    unittest.main()


