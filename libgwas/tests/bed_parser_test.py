#!/usr/bin/env python
import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../../")
    sys.path.insert(0, "../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")
from libgwas import bed_parser
from libgwas.boundary import BoundaryCheck
from libgwas.snp_boundary_check import SnpBoundaryCheck
from libgwas.data_parser import DataParser
from libgwas.pheno_covar import PhenoCovar
import numpy
from libgwas.exceptions import InvalidFrequency
from libgwas.exceptions import TooMuchMissing
from libgwas.exceptions import InvariantVar
import libgwas 

import unittest

from pkg_resources import resource_filename

class TestBase(unittest.TestCase):
    def setUp(self):
        self.missing = "tests/bedfiles/ped_missing"
        self.missing_bed = resource_filename("libgwas", "%s.bed" % (self.missing))
        self.missing_bim = resource_filename("libgwas", "%s.bim" % (self.missing))
        self.missing_fam = resource_filename("libgwas", "%s.fam" % (self.missing))
        self.genotypes = [
            [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
            [1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, 2, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, 1, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, 2, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]
        self.nonmissing = "tests/bedfiles/ped_nomiss"
        self.nonmissing_bed = resource_filename("libgwas", "%s.bed" % (self.nonmissing))
        self.nonmissing_bim = resource_filename("libgwas", "%s.bim" % (self.nonmissing))
        self.nonmissing_fam = resource_filename("libgwas", "%s.fam" % (self.nonmissing))

        self.genotypes_w_missing = [
            [0, 1],
            [1,  1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        self.nonmissing_mapdata = libgwas.get_lines(self.nonmissing_bim, split=True)
        self.missing_mapdata = libgwas.get_lines(self.missing_bim, split=True)

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



class TestBedFilesTPedIndExclusions(TestBase):
    def testCompleteWithExclusions(self):
        DataParser.ind_exclusions = ["1:1", "3:3"]

        genotypes = [
            [1, 0, 1, 0, 0, 1, 0, 0, 1, 0],
            [1, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [2, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [1, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [2, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]

        pc = PhenoCovar()
        ped_parser = bed_parser.Parser(self.nonmissing_fam, self.nonmissing_bim, self.nonmissing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        mapdata = self.nonmissing_mapdata

        index = 0
        for snp in ped_parser:
            for y in pc:
                (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)

                try:
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(genotypes[index], list(genodata.genotypes))
                except TooMuchMissing as e:
                    pass
                except InvalidFrequency as e:
                    pass
            index += 1
        self.assertEqual(7, index)

    def testRegionBoundaryWithExclusions(self):
        DataParser.ind_exclusions = ["1:1", "2:2", "3:3"]

        genotypes = [
            [0, 1, 0, 0, 1, 0, 0, 1, 0],
            [0, 0, 1, 1, 1, 0, 0, 0, 1],
            [1, 0, 0, 0, 2, 1, 1, 0, 0],
            [1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, 0, 0, 1, 2, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]
        BoundaryCheck.chrom = 2
        pc = PhenoCovar()
        ped_parser = bed_parser.Parser(self.nonmissing_fam, self.nonmissing_bim, self.nonmissing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        mapdata = self.nonmissing_mapdata
        index = 4
        for snp in ped_parser:
            for y in pc:
                (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)

                try:
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(genotypes[index], list(genodata.genotypes))

                except TooMuchMissing as e:
                    pass
                except InvalidFrequency as e:
                    pass
            index += 1
        self.assertEqual(7, index)


    def testMissingWithExclusions(self):
        DataParser.ind_exclusions = ["2:2", "3:3"]

        genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, 1],
            [1, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        pc = PhenoCovar()
        ped_parser = bed_parser.Parser(self.missing_fam, self.missing_bim, self.missing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        mapdata = self.missing_mapdata


        self.assertEqual(7, ped_parser.locus_count)
        index = 0
        for snp in ped_parser:
            try:
                for y in pc:
                    (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)

                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(genotypes_w_missing[index], list(genodata.genotypes))
            except TooMuchMissing as e:
                pass
            except InvalidFrequency as e:
                pass
            except InvariantVar as e:
                pass
            index += 1
        self.assertEqual(7, index)

    def testPedWithMissingMxIndExclusionsToo(self):
        pc = PhenoCovar()
        DataParser.ind_exclusions = ["2:2", "3:3"]
        DataParser.ind_miss_tol = 0.5       # We should only lose 1
        ped_parser = bed_parser.Parser(self.missing_fam, self.missing_bim, self.missing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        mapdata = self.missing_mapdata

        genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, 1],
            [1, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        index = 0
        for snp in ped_parser:
            try:
                for y in pc:
                    (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(genotypes_w_missing[index], list(genodata.genotypes))
            except TooMuchMissing as e:
                pass
            except InvalidFrequency as e:
                pass
            except InvariantVar as e:
                pass
            index += 1
        self.assertEqual(7, index)

class TestBedFilesTPed(TestBase):
    # Test to make sure we can load everything
    def testPedComplete(self):
        pc = PhenoCovar()
        ped_parser = bed_parser.Parser(self.nonmissing_fam, self.nonmissing_bim, self.nonmissing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        self.assertEqual(12, ped_parser.ind_count)
        mapdata = self.nonmissing_mapdata

        index = 0
        for snp in ped_parser:
            for y in pc:
                (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)

                try:
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(self.genotypes[index], list(genodata.genotypes))
                except TooMuchMissing as e:
                    pass
                except InvalidFrequency as e:
                    pass
            index += 1
        self.assertEqual(7, index)

    def testTPedPhenoComplete(self):
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        ped_parser = bed_parser.Parser(self.nonmissing_fam, self.nonmissing_bim, self.nonmissing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        self.assertEqual(12, len(pc.covariate_data[0]))
        self.assertEqual(12, len(pc.phenotype_data[0]))
        self.assertEqual(1, len(pc.phenotype_names))
        mapdata = self.nonmissing_mapdata

        index = 0
        for snp in ped_parser:
            for y in pc:
                (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)

                try:
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(self.genotypes[index], list(genodata.genotypes))
                except TooMuchMissing as e:
                    pass
                except InvalidFrequency as e:
                    pass
            index += 1
        self.assertEqual(7, index)


    # Test to make sure we can load everything
    def testPedWithMissingComplete(self):
        pc = PhenoCovar()
        ped_parser = bed_parser.Parser(self.missing_fam, self.missing_bim, self.missing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        mapdata = self.missing_mapdata


        index = 0
        for snp in ped_parser:
            try:
                for y in pc:
                    (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(self.genotypes_w_missing[index], list(genodata.genotypes))
            except TooMuchMissing as e:
                pass
            except InvalidFrequency as e:
                pass
            except InvariantVar as e:
                pass
            index += 1
        self.assertEqual(7, index)

    # Test to make sure we can load everything
    def testPedWithMissingMxIndComplete(self):
        pc = PhenoCovar()
        DataParser.ind_miss_tol = 0.5       # We should only lose 1
        ped_parser = bed_parser.Parser(self.missing_fam, self.missing_bim, self.missing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        mapdata = self.missing_mapdata

        genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1],
            [1,  0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0,  1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0,  2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1,  0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1,  1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0,  1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        index = 0
        for snp in ped_parser:
            try:
                for y in pc:
                    (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(genotypes_w_missing[index], list(genodata.genotypes))
            except TooMuchMissing as e:
                pass
            except InvalidFrequency as e:
                pass
            except InvariantVar as e:
                pass

            index += 1
        self.assertEqual(7, index)

    # Test to make sure we can load everything
    def testPedWithMissingMxSnpComplete(self):
        pc = PhenoCovar()
        DataParser.snp_miss_tol = 0.5       # We should only lose 1
        ped_parser = bed_parser.Parser(self.missing_fam, self.missing_bim, self.missing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        mapdata = self.missing_mapdata
        genotypes_w_missing = [
            [0, 1],
            [1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        index = 0
        for snp in ped_parser:
            try:
                for y in pc:
                    (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(genotypes_w_missing[index], list(genodata.genotypes))
            except TooMuchMissing as e:
                pass
            except InvalidFrequency as e:
                pass
            except InvariantVar as e:
                pass
            index += 1
        self.assertEqual(7, index)

    # Test to make sure we can load everything
    def testPedWithMissingMxBothComplete(self):
        pc = PhenoCovar()
        DataParser.snp_miss_tol = 0.5       # We should only lose 1
        DataParser.ind_miss_tol = 0.5       # We should only lose 1
        ped_parser = bed_parser.Parser(self.missing_fam, self.missing_bim, self.missing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        mapdata = self.missing_mapdata

        genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
            [1,  0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0,  1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0,  2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1,  0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1,  1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0,  1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        index = 0
        valid_loci = 0
        for snp in ped_parser:

            try:
                for y in pc:
                    (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(mapdata[index][0]), snp.chr)
                    self.assertEqual(int(mapdata[index][3]), snp.pos)
                    self.assertEqual(mapdata[index][1], snp.rsid)
                    self.assertEqual(genotypes_w_missing[index], list(genodata.genotypes))
                    valid_loci += 1
            except TooMuchMissing as e:
                pass
            except InvalidFrequency as e:
                pass
            except InvariantVar as e:
                pass
            index += 1
        self.assertEqual(6, valid_loci)
        self.assertEqual(7, index)



    def testPedBoundaryBed(self):
        pc = PhenoCovar()
        DataParser.boundary = BoundaryCheck()
        BoundaryCheck.chrom = 2
        ped_parser = bed_parser.Parser(self.nonmissing_fam, self.nonmissing_bim, self.nonmissing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        pedigree = self.nonmissing_mapdata

        index = 4
        valid_loci = 0
        for snp in ped_parser:
            for y in pc:
                (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)

                try:
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(pedigree[index][0]), snp.chr)
                    self.assertEqual(int(pedigree[index][3]), snp.pos)
                    self.assertEqual(pedigree[index][1], snp.rsid)
                    self.assertEqual(self.genotypes[index], list(genodata.genotypes))
                    valid_loci += 1
                except TooMuchMissing as e:
                    pass
                except InvalidFrequency as e:
                    pass

            index += 1
        self.assertEqual(3, valid_loci)
        self.assertEqual(7, index)

    def testPedSnpBoundaryBed(self):
        pc = PhenoCovar()
        DataParser.boundary = SnpBoundaryCheck(snps=["rs0001-rs0003"])
        BoundaryCheck.chrom = 1
        ped_parser = bed_parser.Parser(self.nonmissing_fam, self.nonmissing_bim, self.nonmissing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        pedigree = self.nonmissing_mapdata

        index = 0
        valid_loci = 0
        self.assertEqual(3, ped_parser.locus_count)
        for snp in ped_parser:
            for y in pc:
                (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)
                try:
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(pedigree[index][0]), snp.chr)
                    self.assertEqual(int(pedigree[index][3]), snp.pos)
                    self.assertEqual(pedigree[index][1], snp.rsid)
                    self.assertEqual(self.genotypes[index], list(genodata.genotypes))
                    valid_loci += 1
                except TooMuchMissing as e:
                    pass
                except InvalidFrequency as e:
                    pass

            index += 1
        self.assertEqual(3, valid_loci)

        # we have selected only the first 3
        self.assertEqual(3, index)

    def testPedSnpBoundary2Bed(self):
        pc = PhenoCovar()
        DataParser.boundary = SnpBoundaryCheck(snps=["rs0005-rs0006"])
        BoundaryCheck.chrom = 2
        ped_parser = bed_parser.Parser(self.nonmissing_fam, self.nonmissing_bim, self.nonmissing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        pedigree = self.nonmissing_mapdata

        index = 4
        self.assertEqual(2, ped_parser.locus_count)
        for snp in ped_parser:
            for y in pc:
                (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)

                try:
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(pedigree[index][0]), snp.chr)
                    self.assertEqual(int(pedigree[index][3]), snp.pos)
                    self.assertEqual(pedigree[index][1], snp.rsid)
                    self.assertEqual(self.genotypes[index], list(genodata.genotypes))
                except TooMuchMissing as e:
                    pass
                except InvalidFrequency as e:
                    pass

            index += 1
        self.assertEqual(6, index)


    def testPedRegionBoundaryTPed(self):
        pc = PhenoCovar()
        DataParser.boundary = SnpBoundaryCheck(snps=["rs0005-rs0006"])
        BoundaryCheck.chrom = 2
        ped_parser = bed_parser.Parser(self.nonmissing_fam, self.nonmissing_bim, self.nonmissing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        pedigree = self.nonmissing_mapdata

        index = 4
        self.assertEqual(2, ped_parser.locus_count)
        for snp in ped_parser:
            for y in pc:
                (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)

                try:
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(pedigree[index][0]), snp.chr)
                    self.assertEqual(int(pedigree[index][3]), snp.pos)
                    self.assertEqual(pedigree[index][1], snp.rsid)
                    self.assertEqual(self.genotypes[index], list(genodata.genotypes))
                except TooMuchMissing as e:
                    pass
                except InvalidFrequency as e:
                    pass

            index += 1
        self.assertEqual(6, index)

    def testPedRegionBoundaryWithExclusionsTPed(self):
        pc = PhenoCovar()
        DataParser.boundary = SnpBoundaryCheck(snps=["rs0005-rs0007"])
        DataParser.boundary.LoadExclusions(snps=["rs0007"])
        BoundaryCheck.chrom = 2
        ped_parser = bed_parser.Parser(self.nonmissing_fam, self.nonmissing_bim, self.nonmissing_bed)
        ped_parser.load_fam(pc)
        ped_parser.load_bim(map3=False)
        ped_parser.load_genotypes()

        pedigree = self.nonmissing_mapdata

        index = 4
        self.assertEqual(2, ped_parser.locus_count)
        for snp in ped_parser:
            for y in pc:
                (pheno, covars, nonmissing) = y.get_variables(snp.missing_genotypes)

                try:
                    genodata = snp.get_genotype_data(nonmissing)
                    self.assertEqual(int(pedigree[index][0]), snp.chr)
                    self.assertEqual(int(pedigree[index][3]), snp.pos)
                    self.assertEqual(pedigree[index][1], snp.rsid)
                    self.assertEqual(self.genotypes[index], list(genodata.genotypes))
                except TooMuchMissing as e:
                    pass
                except InvalidFrequency as e:
                    pass
            index += 1
        self.assertEqual(6, index)


if __name__ == "__main__":
    unittest.main()
