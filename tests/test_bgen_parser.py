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
from libgwas.pheno_covar import PhenoCovar, PhenoIdFormat
from libgwas import bgen_parser
from libgwas.boundary import BoundaryCheck
from pkg_resources import resource_filename
import libgwas.impute_parser

# import gzip

base_freq = [0.99, 0.75, 0.7,0.8, 0.65, 0.7, 0.85, 0.7, 0.7, 0.3]


class TestBase(unittest.TestCase):
    def setUp(self):
        self.allele_1 = list("AAACCCGGGTCGTGTATACC")
        self.allele_2 = list("CGTGTATACCAAACCCGGGT")

        self.nomissing = resource_filename("tests", 'bedfiles/test.bgen')
        self.nomissing_sample = resource_filename("tests", 'bedfiles/test.bgen.sample')


        self.additive_encoding = [
            [0.000000, 0.127169, 0.000000, 0.026108, 0.151843, 0.104356, 0.040452, 0.106310, 0.026459, 0.403708,
            0.000000, 0.000000],
            [0.379370, 0.568566, 0.938415, 0.404898, 0.627787, 0.835401, 0.230304, 0.406455, 0.384497, 0.456489,
            0.616754, 0.345785],
            [0.819898, 1.012375, 0.567895, 0.703471, 0.915068, 0.590921, 0.432792, 0.760784, 0.630350, 0.562997,
            0.971618, 0.248157],
            [0.607736, 0.702739, 0.381674, 0.360235, 0.203494, 0.641962, 0.351629, 0.203464, 0.159060, 0.550835,
            0.597894, 0.145480],
            [0.580652, 0.666178, 0.798642, 0.815656, 0.881361, 0.369909, 1.057939, 0.671199, 0.980896, 0.551324,
            0.655192, 0.928801],
            [0.677592, 0.623163, 0.691890, 0.922942, 0.992157, 0.338140, 0.410758, 0.317647, 0.689769, 0.341741,
            0.765316, 0.213367],
            [0.402716, 0.352697, 0.111910, 0.208103, 0.022477, 0.166033, 0.568475, 0.125948, 0.117266, 0.557259,
            0.403021, 0.228489],
            [0.332448, 0.555749, 0.589624, 0.275273, 0.762158, 0.189563, 0.594308, 0.625116, 0.263355, 0.747143,
            0.453483, 0.549889],
            [0.659739, 0.697078, 0.876234, 0.766857, 0.345663, 0.712200, 0.643534, 0.642817, 0.706020, 0.452323,
            0.394064, 0.606271],
            [1.424872, 1.197848, 1.489326, 1.362325, 1.172427, 1.503670, 1.514473, 1.507042, 1.345663, 1.486091,
            1.292332, 1.445747],
            [0.284001, 0.000000, 0.000000, 0.208789, 0.000000, 0.068513, 0.075715, 0.050599, 0.014786, 0.000000,
            0.000000, 0.000000],
            [0.506264, 0.224201, 0.860700, 0.196689, 0.476310, 0.524498, 0.719203, 0.851087, 0.441306, 0.518334,
            0.532097, 0.626398],
            [0.774090, 0.838102, 0.662806, 0.449104, 0.685237, 0.807706, 0.912200, 0.833310, 0.515007, 0.628107,
            0.591669, 0.638407],
            [0.411810, 0.125216, 0.350500, 0.838285, 0.281407, 0.657359, 0.470710, 0.187610, 0.432212, 0.587961,
            0.171679, 0.731609],
            [0.513161, 0.909300, 0.455039, 1.007599, 1.106294, 0.633326, 0.528801, 0.509804, 0.962905, 0.570611,
            0.773098, 0.491997],
            [0.763394, 0.872496, 0.506065, 0.640696, 0.291096, 0.655985, 0.750210, 0.628885, 0.465293, 0.513451,
            0.625498, 0.737087],
            [0.345403, 0.168490, 0.484215, 0.169299, 0.511498, 0.295995, 0.063813, 0.405005, 0.810697, 0.227893,
            0.079698, 0.400641],
            [0.296605, 0.445899, 0.180896, 0.886610, 0.728588, 0.685908, 0.701930, 0.495094, 0.417212, 0.123598,
            0.828504, 0.589837],
            [0.297200, 0.624506, 0.787900, 0.609308, 0.674388, 0.695109, 0.691096, 0.537652, 0.870207, 0.834089,
            0.492699, 0.322377],
            [1.807797, 1.595346, 1.613611, 1.205295, 1.541299, 1.484199, 1.312688, 1.560998, 1.517891, 1.587198,
            1.231907, 1.528298]
        ]
        self.phenotypes = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5, 0.9, 1.0]
        self.sex = [1,2,1,1,2,2,1,1,2,2,2,1]

        self.positions = [13970,15367,16764,18161,19558,20955,22352,23749,25146,26543,27940,29337,30734,32131,33528,34925,36322,37719,39116,40513]
        self.mafs = [
            0.0411001754788,
            0.258113349101,
            0.342346964726,
            0.20442511635,
            0.373239490349,
            0.291020065614,
            0.136016378017,
            0.247421225299,
            0.312616667938,
            0.697575722896,
            0.0292668039979,
            0.269878690776,
            0.34732267745,
            0.218598204522,
            0.352580682078,
            0.310423183541,
            0.165110246433,
            0.265861753262,
            0.309855420768,
            0.749438595153
        ]
        self.rsids = "rs1320,rs13267,rs132134,rs132201,rs132268,rs132335,rs132402,rs132469,rs132536,rs132603,rs132670,rs132737,rs132804,rs132871,rs132938,rs1321005,rs1321072,rs1321139,rs1321206,rs1321273".split(",")
        self.ind_ids = 'ID0001,ID20002,ID0003,ID0004,IID0005,0006,ID0007,ID0008,ID0009,ID0010,ID0011,IID0012'.split(',')
        # For now, we aren't calculating the info scores and it always returns 1.0
        self.info = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]


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
        self.encoding = bgen_parser.encoding
        self.parser_info_thresh = bgen_parser.Parser.info_threshold
        bgen_parser.Parser.info_threshold = 0.0

        self.raw = numpy.zeros((20, 12, 3))

        self.orig_id_encoding = PhenoCovar.id_encoding

        # This is required for the current bgen sample file format
        PhenoCovar.id_encoding = PhenoIdFormat.FID

    def tearDown(self):
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
        bgen_parser.encoding = self.encoding
        bgen_parser.Parser.info_threshold = self.parser_info_thresh
        PhenoCovar.id_encoding = self.orig_id_encoding

class TestBGenBasics(TestBase):
    def testWithSample(self):

        pc = PhenoCovar()
        parser = bgen_parser.Parser(self.nomissing, self.nomissing_sample)

        parser.load_family_details(pc)
        parser.load_genotypes()
        pc.load_phenofile(open(self.nomissing_sample), names=['plink_pheno'], sample_file=True)
        pc.load_covarfile(open(self.nomissing_sample), names=['sex'], sample_file=True)

        idx = 0
        for id in self.ind_ids:
            self.assertTrue(id in pc.pedigree_data)
            self.assertEqual(self.phenotypes[idx], pc.phenotype_data[0][idx])
            self.assertEqual(self.sex[idx], pc.covariate_data[0][idx])
            idx += 1
        self.assertEqual(12, idx)

    def testWithoutSample(self):
        pc = PhenoCovar()
        parser = bgen_parser.Parser(self.nomissing)

        parser.load_family_details(pc)
        parser.load_genotypes()
        idx = 0
        for id in ["ID%s" % str(x).zfill(4) for x in range(1,12)]:
            self.assertTrue(id in pc.pedigree_data)

    def testAdditiveGeno(self):
        pc = PhenoCovar()
        parser = bgen_parser.Parser(self.nomissing)

        parser.load_family_details(pc)
        parser.load_genotypes()
        idx = 0
        for snp in parser:
            self.assertEqual(self.positions[idx], snp.pos)
            for ind in range(snp.genotype_data.shape[0]):
                self.assertAlmostEqual(self.additive_encoding[idx][ind], snp.genotype_data[ind], places=4)
            idx += 1

if __name__ == "__main__":
    unittest.main()

