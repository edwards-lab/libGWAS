#!/usr/bin/env python
import sys

if "DEBUG" in sys.argv:
    sys.path.insert(0, "../../")
    sys.path.insert(0, "../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import unittest
import os
from libgwas.data_parser import DataParser
from libgwas.pheno_covar import PhenoCovar
from libgwas.vcf_parser import Parser
from libgwas.boundary import BoundaryCheck
from libgwas.snp_boundary_check import SnpBoundaryCheck
from pkg_resources import resource_filename

class TestBase(unittest.TestCase):
    def setUp(self):
        self.nonmissing = resource_filename("tests", "bedfiles/nomiss.vcf")
        self.nonmissinggz = resource_filename("tests", "bedfiles/nomiss.vcf.gz")
        self.genotypes = [
            [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
            [1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, 2, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, 1, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, 2, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]
        self.missing = resource_filename("tests", "bedfiles/miss.vcf")
        self.missinggz = resource_filename("tests", "bedfiles/miss.vcf.gz")

        self.genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1],
            [1,  1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, -1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, -1, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, -1, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, -1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, -1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]

        self.nonmissing_mapdata = []
        for line in open(self.nonmissing):
            if line[0] != "#":
                words = line.strip().split()
                # CHR   Pos ID  REF ALT
                self.nonmissing_mapdata.append(words[0:5])
        self.missing_mapdata = []
        for line in open(self.missing):
            if line[0] != "#":
                words = line.strip().split()
                # CHR   Pos ID  REF ALT
                self.missing_mapdata.append(words[0:5])
        self.phenotypes     = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5, 0.9, 1.0]
        self.sex            = [1,1,2,2,1,1,1,1,2,2,1,1]

        self.chrom          = BoundaryCheck.chrom
        self.boundary       = DataParser.boundary
        self.min_maf        = DataParser.min_maf
        self.max_maf        = DataParser.max_maf
        self.snp_miss_tol   = DataParser.snp_miss_tol
        self.ind_miss_tol   = DataParser.ind_miss_tol
        self.sex_as_covar = PhenoCovar.sex_as_covariate
        self.has_sex        = DataParser.has_sex
        self.has_pheno      = DataParser.has_pheno
        self.has_parents    = DataParser.has_parents
        self.has_fid        = DataParser.has_fid
        self.has_liability  = DataParser.has_liability

        DataParser.boundary = BoundaryCheck()

    def tearDown(self):
        PhenoCovar.sex_as_covariate = self.sex_as_covar
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

class TestVcfFiles(TestBase):
    def testNonmissingBasics(self):
        pc = PhenoCovar()
        parser = Parser(self.nonmissing, data_field='GT')
        parser.init_subjects(pc)
        parser.load_genotypes()

        self.assertEqual(0, len(pc.covariate_labels))
        self.assertEqual(0, len(pc.phenotype_names))
        mapdata = self.nonmissing_mapdata

        index = 0
        for snp in parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][1]), snp.pos)
            self.assertEqual(mapdata[index][2], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testGzBasics(self):
        pc = PhenoCovar()
        parser = Parser(self.nonmissinggz, data_field='GT')
        parser.init_subjects(pc)
        parser.load_genotypes()

        self.assertEqual(0, len(pc.covariate_labels))
        self.assertEqual(0, len(pc.phenotype_names))
        mapdata = self.nonmissing_mapdata

        index = 0
        for snp in parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][1]), snp.pos)
            self.assertEqual(mapdata[index][2], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testMissingBasics(self):
        pc = PhenoCovar()
        parser = Parser(self.missing, data_field='GT')
        parser.init_subjects(pc)
        parser.load_genotypes()

        self.assertEqual(0, len(pc.covariate_labels))
        self.assertEqual(0, len(pc.phenotype_names))
        mapdata = self.missing_mapdata

        index = 0
        for snp in parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][1]), snp.pos)
            self.assertEqual(mapdata[index][2], snp.rsid)
            self.assertEqual(self.genotypes_w_missing[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testMissingMxIndComplete(self):
        pc = PhenoCovar()
        DataParser.ind_miss_tol = 0.5       # We should only lose 1
        parser = Parser(self.missing, data_field='GT')
        parser.init_subjects(pc)
        parser.load_genotypes()

        mapdata = self.missing_mapdata

        genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1],
            [1,  0,  0,  0,  1,  1,  1,  0,  0,  0, 1],
            [0,  1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0,  2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1,  0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1,  1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0,  1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        index = 0
        for snp in parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][1]), snp.pos)
            self.assertEqual(mapdata[index][2], snp.rsid)
            self.assertEqual(genotypes_w_missing[index],
                             list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testMissingMxSnpComplete(self):
        pc = PhenoCovar()
        DataParser.snp_miss_tol = 0.5       # We should only lose 1
        parser = Parser(self.missing, data_field='GT')
        parser.init_subjects(pc)
        parser.load_genotypes()

        mapdata = self.missing_mapdata

        index = 1
        for snp in parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][1]), snp.pos)
            self.assertEqual(mapdata[index][2], snp.rsid)
            self.assertEqual(self.genotypes_w_missing[index],
                             list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    # Test to make sure we can load everything
    def testMissingMxBothComplete(self):
        pc = PhenoCovar()
        DataParser.snp_miss_tol = 0.5       # We should only lose 1
        DataParser.ind_miss_tol = 0.5       # We should only lose 1
        parser = Parser(self.missing, data_field='GT')
        parser.init_subjects(pc)
        parser.load_genotypes()

        mapdata = self.missing_mapdata

        index = 1
        genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
            [1,  0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0,  1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0,  2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1,  0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1,  1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0,  1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        for snp in parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][1]), snp.pos)
            self.assertEqual(mapdata[index][2], snp.rsid)
            self.assertEqual(genotypes_w_missing[index],
                             list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

if __name__ == "__main__":
    unittest.main()
