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
from libgwas import impute_parser
from libgwas.boundary import BoundaryCheck


import gzip

base_freq = [0.99, 0.75, 0.7,0.8, 0.65, 0.7, 0.85, 0.7, 0.7, 0.3]


class TestBase(unittest.TestCase):
    def setUp(self):
        self.allele_1 = list("AAACCCGGGTCGTGTATACC")
        self.allele_2 = list("CGTGTATACCAAACCCGGGT")

        self.WriteTestFiles()

        self.phenotypes = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5, 0.9, 1.0]
        self.sex = [1,2,1,1,2,2,1,1,2,2,2,1]

        self.gen_ext = impute_parser.Parser.gen_ext
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
        self.encoding = impute_parser.encoding
        self.compression = DataParser.compressed_pedigree
        DataParser.compressed_pedigree = True
        self.parser_info_thresh = impute_parser.Parser.info_threshold
        impute_parser.Parser.info_threshold = 0.0

    def tearDown(self):
        os.remove(self.fam_file)
        os.remove(self.gen_file)
        os.remove(self.gen_file2)
        os.remove(self.uncmp_1)
        os.remove(self.uncmp_2)
        os.remove(self.info_file1)
        os.remove(self.info_file2)

        impute_parser.Parser.gen_ext = self.gen_ext
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
        impute_parser.encoding = self.encoding
        DataParser.compressed_pedigree = self.compression
        impute_parser.Parser.info_threshold = self.parser_info_thresh

    def WriteTestFiles(self, prefix = "__test_imputed"):

        self.fam_file = "%s.gen_samples" % (prefix)
        fam_file = open(self.fam_file, 'w')
        print >> fam_file, """ID_1 ID_2 missing father mother sex plink_pheno
    0 0 0 D D D B
    ID0001 FAM001 0 0 0 1 0.1
    ID0002 FAM002 0 0 0 2 0.4
    ID0003 FAM003 0 0 0 1 1.0
    ID0004 FAM004 0 0 0 1 0.5
    ID0005 FAM005 0 0 0 2 0.9
    ID0006 FAM006 0 0 0 2 1.0
    ID0007 FAM007 0 0 0 1 0.1
    ID0008 FAM008 0 0 0 1 0.4
    ID0009 FAM009 0 0 0 2 1.0
    ID0010 FAM010 0 0 0 2 0.5
    ID0011 FAM011 0 0 0 2 0.9
    ID0012 FAM012 0 0 0 1 1.0"""
        fam_file.close()
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

        self.gen_file = "%s.gen.gz" % (prefix)
        self.gen_file2 = "%s-2.gen.gz" % (prefix)
        self.info_file1 = "%s.info" % (prefix)
        self.info_file2 = "%s-2.info" % (prefix)
        self.uncmp_1 = "%s.gen" % (prefix)
        self.uncmp_2 = "%s-2.gen" % (prefix)
        gen_file = gzip.open(self.gen_file, 'wb')
        uncmp_file = open(self.uncmp_1, 'w')
        idx = 0
        self.additive_encoding = numpy.zeros((20, 12))
        self.dominant_encoding = numpy.zeros((20, 12))
        self.recessive_encoding = numpy.zeros((20, 12))
        self.raw = numpy.zeros((20, 12, 3))
        self.positions = []
        self.mafs = []
        self.rsids = []

        info = [0.1 * (x%10) + 0.05 for x in range(0, 20)]
        certainty = [0.2 * (x%5) + 0.05 for x in range(0, 20)]

        info_file = open(self.info_file1, 'w')
        print >> info_file, "snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0"

        for base in base_freq:
            f = numpy.random.normal(loc=base, scale=0.1, size=12)
            f[f>1.0] = 1.0
            f[f<0] = 0.0
            maf = 1.0 - f
            AA = f * f
            Aa = 2 * f * maf
            aa = maf * maf

            self.raw[idx] = numpy.hstack((AA.reshape(-1,1), Aa.reshape(-1, 1), aa.reshape(-1, 1)))
            self.additive_encoding[idx] = Aa + 2*aa
            self.mafs.append(numpy.mean(self.additive_encoding[idx]/2))
            self.dominant_encoding[idx] = Aa + aa
            self.recessive_encoding[idx] = aa
            line = numpy.hstack((AA.reshape(-1, 1), Aa.reshape(-1, 1), aa.reshape(-1, 1)))
            self.positions.append((10+idx) * 1397)
            self.rsids.append("rs132%d" % (idx * 67))
            print >> gen_file, "\t".join([
                                            "--",
                                            self.rsids[-1],
                                            str(self.positions[-1]),
                                            self.allele_1[idx],
                                            self.allele_2[idx]]) + \
                               "\t" + \
                               "\t".join(["%0.6f" % (x) for x in line.reshape(-1)])
            print >> uncmp_file, "\t".join([
                                            "--",
                                            self.rsids[-1],
                                            str(self.positions[-1]),
                                            self.allele_1[idx],
                                            self.allele_2[idx]]) + \
                               "\t" + \
                               "\t".join(["%0.6f" % (x) for x in line.reshape(-1)])
            print >> info_file, " ".join([
                                            "--",
                                            self.rsids[-1],
                                            str(self.positions[-1]),
                                            str(numpy.mean(self.additive_encoding[idx]/2)),
                                            str(info[idx]),
                                            str(certainty[idx]),
                                            "0",
                                            "-1",
                                            "-1",
                                            "-1"])


            idx += 1
        gen_file.close()
        uncmp_file.close()
        info_file.close()
        gen_file = gzip.open(self.gen_file2, 'wb')
        uncmp_file = open(self.uncmp_2, 'w')
        info_file = open(self.info_file2, 'w')
        print >> info_file, "snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0"
        for base in base_freq:
            f = numpy.random.normal(loc=base, scale=0.1, size=12)
            f[f>1.0] = 1.0
            f[f<0] = 0.0
            maf = 1.0 - f
            AA = f * f
            Aa = 2 * f * maf
            aa = maf * maf
            self.raw[idx] = numpy.hstack((AA.reshape(-1,1), Aa.reshape(-1, 1), aa.reshape(-1, 1)))
            self.additive_encoding[idx] = Aa + 2*aa
            self.mafs.append(numpy.mean(self.additive_encoding[idx]/2))
            self.dominant_encoding[idx] = Aa + aa
            self.recessive_encoding[idx] = aa
            line = numpy.hstack((AA.reshape(-1, 1), Aa.reshape(-1, 1), aa.reshape(-1, 1)))
            self.positions.append((10+idx) * 1397)
            self.rsids.append("rs132%d" % (idx * 67))
            print >> gen_file, "\t".join([
                                            "--",
                                            self.rsids[-1],
                                            str(self.positions[-1]),
                                            self.allele_1[idx],
                                            self.allele_2[idx]]) + \
                               "\t" + \
                               "\t".join(["%0.4f" % (x) for x in line.reshape(-1)])
            print >> uncmp_file, "\t".join([
                                            "--",
                                            self.rsids[-1],
                                            str(self.positions[-1]),
                                            self.allele_1[idx],
                                            self.allele_2[idx]]) + \
                               "\t" + \
                               "\t".join(["%0.4f" % (x) for x in line.reshape(-1)])
            print >> info_file, " ".join([
                                            "--",
                                            self.rsids[-1],
                                            str(self.positions[-1]),
                                            str(numpy.mean(self.additive_encoding[idx]/2)),
                                            str(info[idx]),
                                            str(certainty[idx]),
                                            "0",
                                            "-1",
                                            "-1",
                                            "-1"])
            idx += 1
        self.impute_parser = impute_parser.Parser(self.fam_file, [self.gen_file], chroms = ["3"])

class TestImputedBasics(TestBase):
    def testFamilyData(self):
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0
        for id in self.ind_ids:
            self.assertTrue(id in pc.pedigree_data)
            self.assertEqual(self.phenotypes[idx], pc.phenotype_data[0][idx])
            self.assertEqual(self.sex[idx], pc.covariate_data[0][idx])
            idx += 1

    def testChromosomes(self):
        impute_parser.encoding = impute_parser.Encoding.Raw
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"])
        chromosomes = ['3']*10+['4']*10
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(chromosomes[idx], snp.chr)

            idx += 1
        self.assertEqual(20, idx)
    def testRawValuesUncompressed(self):
        DataParser.compressed_pedigree = False
        impute_parser.Parser.gen_ext = "gen"
        impute_parser.encoding = impute_parser.Encoding.Raw
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.uncmp_1, self.uncmp_2], chroms = ["3", "4"])
        chromosomes = (['3']*10+['4']*10)*2
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(chromosomes[idx], snp.chr)

            for i in range(0, len(self.raw[idx])):
                for j in [0, 1, 2]:
                    self.assertAlmostEqual(self.raw[idx][i][j], snp.genotype_data[i][j], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testAdditiveValuesUncompressed(self):
        impute_parser.Parser.gen_ext = "gen"
        DataParser.compressed_pedigree = False
        impute_parser.encoding = impute_parser.Encoding.Additive
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.uncmp_1, self.uncmp_2], chroms = ["3", "4"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.additive_encoding[idx])):
                self.assertAlmostEqual(self.additive_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testRawValues(self):
        impute_parser.encoding = impute_parser.Encoding.Raw
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"])
        chromosomes = (['3']*10+['4']*10)*2
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(chromosomes[idx], snp.chr)

            for i in range(0, len(self.raw[idx])):
                for j in [0, 1, 2]:
                    self.assertAlmostEqual(self.raw[idx][i][j], snp.genotype_data[i][j], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testAdditiveValues(self):
        impute_parser.encoding = impute_parser.Encoding.Additive
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.additive_encoding[idx])):
                self.assertAlmostEqual(self.additive_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testAdditiveValuesFilterInfo(self):
        impute_parser.Parser.info_threshold = 0.4
        impute_parser.encoding = impute_parser.Encoding.Additive
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file], chroms = ["3"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 4

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.additive_encoding[idx])):
                self.assertAlmostEqual(self.additive_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(10, idx)

    def testDominantValues(self):
        impute_parser.encoding = impute_parser.Encoding.Dominant
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.dominant_encoding[idx])):
                self.assertAlmostEqual(self.dominant_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testInfoFileUse(self):
        # We'll give it an invalid gen_ext so that we can be certain that it's using the files provided
        impute_parser.Parser.gen_ext='asdf'
        impute_parser.encoding = impute_parser.Encoding.Dominant
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"], info_files=[self.info_file1, self.info_file2])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.dominant_encoding[idx])):
                self.assertAlmostEqual(self.dominant_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testRecessiveValues(self):
        impute_parser.encoding = impute_parser.Encoding.Recessive
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.recessive_encoding[idx])):
                self.assertAlmostEqual(self.recessive_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testGenotypeValues(self):
        impute_parser.encoding = impute_parser.Encoding.Genotype
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.raw[idx])):
                genotype = 2
                AA, Aa, aa = self.raw[idx][i]
                if Aa >= AA and Aa >= aa:
                    genotype = 1
                elif AA >= Aa and AA >= aa:
                    genotype = 0
                self.assertAlmostEqual(genotype, snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testMAF(self):
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file], chroms = ["3"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            maf = numpy.mean(self.mafs[idx])
            if maf > 0.5:
                maf = 1.0 - maf
            self.assertAlmostEqual(maf, snp.maf, places=3)
            idx += 1
        self.assertEqual(10, idx)

    def testAlelles(self):
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file], chroms = ["3"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(snp.major_allele, self.allele_1[idx])
            self.assertEqual(snp.minor_allele, self.allele_2[idx])
            self.assertEqual(snp.rsid, self.rsids[idx])
            idx += 1
        self.assertEqual(10, idx)

    def testLongerList(self):
        impute_parser.encoding = impute_parser.Encoding.Recessive
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2]*3, chroms=["3", "4"]*3)
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        values = numpy.vstack((self.recessive_encoding, self.recessive_encoding, self.recessive_encoding))
        positions = self.positions * 3

        for snp in parser:
            self.assertEqual(positions[idx], snp.pos)
            for i in range(0, len(values[idx])):
                self.assertAlmostEqual(values[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(60, idx)

    def testFilteredInd(self):
        DataParser.ind_exclusions = self.ind_ids[0:2]
        impute_parser.encoding = impute_parser.Encoding.Recessive
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(2, len(self.recessive_encoding[idx])):
                self.assertAlmostEqual(self.recessive_encoding[idx][i], snp.genotype_data[i-2], places=3)
            idx += 1
        self.assertEqual(20, idx)

    def testBoundaried(self):
        BoundaryCheck.chrom = 3
        DataParser.boundary = BoundaryCheck(bp=[0, 20900])
        impute_parser.encoding = impute_parser.Encoding.Recessive
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = [3, 4])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.recessive_encoding[idx])):
                self.assertAlmostEqual(self.recessive_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(5, idx)

    def testBoundariedUpper(self):
        BoundaryCheck.chrom = 3
        DataParser.boundary = BoundaryCheck(bp=[21000, 50000])
        impute_parser.encoding = impute_parser.Encoding.Recessive
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = [3, 4])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 6

        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.recessive_encoding[idx])):
                self.assertAlmostEqual(self.recessive_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1

    def testBoundariedMiddle(self):
        BoundaryCheck.chrom = 4
        DataParser.boundary = BoundaryCheck(bp=[30734, 33528])
        impute_parser.encoding = impute_parser.Encoding.Recessive
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = [3, 4])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0
        dropped = 0
        for snp in parser:
            while self.positions[idx] < 30734 or self.positions[idx] > 33528:
                idx += 1
                dropped += 1
            self.assertEqual(self.positions[idx], snp.pos)
            for i in range(0, len(self.recessive_encoding[idx])):
                self.assertAlmostEqual(self.recessive_encoding[idx][i], snp.genotype_data[i], places=3)
            idx += 1
        self.assertEqual(12, dropped)


    def testFilterMAF(self):
        DataParser.min_maf = 0.45
        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file], chroms = ["3"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0
        for snp in parser:
            while numpy.mean(self.mafs[idx]) < DataParser.min_maf:
                idx += 1
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(snp.major_allele, self.allele_1[idx])
            self.assertEqual(snp.minor_allele, self.allele_2[idx])
            self.assertEqual(snp.rsid, self.rsids[idx])
            idx += 1

    def testFilterSNP(self):
        DataParser.boundary.LoadExclusions(snps=["rs132670", "rs132938"])

        pc = PhenoCovar()
        parser = impute_parser.Parser(self.fam_file, [self.gen_file, self.gen_file2], chroms = ["3", "4"])
        parser.load_family_details(pc)
        parser.load_genotypes()

        idx = 0
        dropped = 0
        for snp in parser:
            while self.rsids[idx] in DataParser.boundary.ignored_rs:
                dropped += 1
                idx += 1
            self.assertEqual(self.positions[idx], snp.pos)
            self.assertEqual(snp.major_allele, self.allele_1[idx])
            self.assertEqual(snp.minor_allele, self.allele_2[idx])
            self.assertEqual(snp.rsid, self.rsids[idx])
            idx += 1
        self.assertEqual(2, dropped)





if __name__ == "__main__":
    unittest.main()


