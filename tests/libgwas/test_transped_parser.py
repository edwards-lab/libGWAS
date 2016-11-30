#!/usr/bin/env python
import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../../")
    sys.path.insert(0, "../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import unittest
import os
from libgwas.data_parser import DataParser
from libgwas.pheno_covar import PhenoCovar
from libgwas.transposed_pedigree_parser import Parser as TransposedPedigreeParser
from libgwas.boundary import BoundaryCheck
from libgwas.snp_boundary_check import SnpBoundaryCheck

class TestBase(unittest.TestCase):
    def setUp(self):
        self.WriteTestFiles()

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

        for file in self.filenames:
            os.remove(file)
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



    def WriteTestFiles(self, prefix = "__test_pedigree"):
        self.filenames = []
        self.tfam_filename = "%s.tfam" % (prefix)
        self.filenames.append(self.tfam_filename)

        f = open(self.tfam_filename, "w")
        f.write("""1 1 0 0 1 0.1
2 2 0 0 1 0.4
3 3 0 0 2 1.0
4 4 0 0 2 0.5
5 5 0 0 1 0.9
6 6 0 0 1 1.0
7 7 0 0 1 0.1
8 8 0 0 1 0.4
9 9 0 0 2 1.0
10 10 0 0 2 0.5
11 11 0 0 1 0.9
12 12 0 0 1 1.0

""")


        f.close()

        self.tped_filename = "%s.tped" % (prefix)
        self.filenames.append(self.tped_filename)
        f = open(self.tped_filename, "w")
        f.write("""1 rs0001 0 500 A A A C A A A A A C A A A A A C A A A A A C A A
1 rs0002 0 10000 G T G T G G G G G G G T G T G T G G G G G G G T
1 rs0003 0 25000 A A G G A G A G A A A A A A G G A G A G A A A A
1 rs0004 0 45000 G G C G C C C G C G G G G G C G C C C G C G G G
2 rs0005 0 750 C T T T C C C T C C C C C T T T C C C T C C C C
2 rs0006 0 10000 G T G G G T G G G G G G G G G G G T G G G G G G
2 rs0007 0 25000 T T C T C T T T T T T T T T C T C T T T T T T T
""")
        f.close()
        self.tped1_alleles = [["A","C"],["G","T"], ["A","G"],["G","C"], ["C","T"],["G","T"],["T","C"]]
        self.hetero_freq_tped = [0.3333, 0.5, 0.3333, 0.41667, 0.3333, 0.25,0.3333]
        self.miss_tped_filename = "%s-missing.tped" % (prefix)
        self.filenames.append(self.miss_tped_filename)
        f = open(self.miss_tped_filename, "w")
        f.write("""1 rs0001 0 500 A A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 A C
1 rs0002 0 10000 G T G T G G G G G G G T G T G T G G G G G G G T
1 rs0003 0 25000 A A 0 0 A G A G A A A A A A G G A G A G A A A A
1 rs0004 0 45000 G G 0 0 C C C G C G G G G G C G C C C G C G G G
2 rs0005 0 750 C T 0 0 C C C T C C C C C T T T C C C T C C C C
2 rs0006 0 10000 G T 0 0 G T G G G G G G G G G G G T G G G G G G
2 rs0007 0 25000 T T 0 0 C T T T T T T T T T C T C T T T T T T T
""")
        self.hetero_freq_misstped = [0.0833, 0.5, 0.3333, 0.41667, 0.3333, 0.25,0.3333]

        self.genotypes = [
            [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
            [1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, 2, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, 1, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, 2, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]
        self.genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1],
            [1,  1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, -1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, -1, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, -1, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, -1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, -1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        f.close()
        self.misssnp_tped_filename = "%s-miss_snps.tped" % (prefix)
        self.filenames.append(self.misssnp_tped_filename)
        f = open(self.misssnp_tped_filename, "w")
        f.write("""1 rs0001 0 -500 A A A C A A A A A C A A A A A C A A A A A C A A
1 rs0002 0 -10000 G T G T G G G G G G G T G T G T G G G G G G G T
1 rs0003 0 25000 A A G G A G A G A A A A A A G G A G A G A A A A
1 rs0004 0 45000 G G C G C C C G C G G G G G C G C C C G C G G G
2 rs0005 0 750 C T T T C C C T C C C C C T T T C C C T C C C C
2 rs0006 0 10000 G T G G G T G G G G G G G G G G G T G G G G G G
2 rs0007 0 25000 T T C T C T T T T T T T T T C T C T T T T T T T
""")
        f.close()
class TestPedFilesTPedIndExclusions(TestBase):
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
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(genotypes[index], list(snp.genotype_data))
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
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]
        index = 4
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(genotypes[index], list(snp.genotype_data))
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
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.miss_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.miss_tped_filename).readlines()]


        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(genotypes_w_missing[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testPedWithMissingMxIndExclusionsToo(self):
        pc = PhenoCovar()
        DataParser.ind_exclusions = ["2:2", "3:3"]
        DataParser.ind_miss_tol = 0.5       # We should only lose 1
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.miss_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.miss_tped_filename).readlines()]

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
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(genotypes_w_missing[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

class TestPedFilesTPed(TestBase):
    # Test to make sure we can load everything
    def testPedComplete(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()


        self.assertEqual(12, ped_parser.ind_count)
        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testTPedPhenoComplete(self):
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)

        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        self.assertEqual(12, len(pc.covariate_data[0]))
        self.assertEqual(12, len(pc.phenotype_data[0]))
        self.assertEqual(1, len(pc.phenotype_names))
        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        self.assertEqual(7, ped_parser.locus_count)
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)


    # Test to make sure we can load everything
    def testPedWithMissingComplete(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.miss_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.miss_tped_filename).readlines()]

        self.assertEqual(7, ped_parser.locus_count)

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes_w_missing[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    # Test to make sure we can load everything
    def testPedWithMissingMxIndComplete(self):
        pc = PhenoCovar()
        DataParser.ind_miss_tol = 0.5       # We should only lose 1
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.miss_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.miss_tped_filename).readlines()]
        genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1],
            [1,  0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0,  1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0,  2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1,  0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1,  1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0,  1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        self.assertEqual(7, ped_parser.locus_count)

        index = 0
        loci = ped_parser.get_loci()
        for snp in loci:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            index += 1

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(genotypes_w_missing[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    # Test to make sure we can load everything
    def testPedWithMissingMxSnpComplete(self):
        pc = PhenoCovar()
        DataParser.snp_miss_tol = 0.5       # We should only lose 1
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.miss_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.miss_tped_filename).readlines()]

        genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
            [1,  1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, -1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, -1, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, -1, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, -1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, -1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]

        hetero_freq_tped = [0.3636, 0.5, 0.3636, 0.4545, 0.3636, 0.2727, 0.2727]

        self.assertEqual(6, ped_parser.locus_count)
        index = 1
        loci = ped_parser.get_loci()
        for snp in loci:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertAlmostEqual(hetero_freq_tped[index], snp.hetero_freq, places=4)
            index += 1

        index = 1
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(genotypes_w_missing[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    # Test to make sure we can load everything
    def testPedWithMissingMxBothComplete(self):
        pc = PhenoCovar()
        DataParser.snp_miss_tol = 0.5       # We should only lose 1
        DataParser.ind_miss_tol = 0.5       # We should only lose 1
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.miss_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.miss_tped_filename).readlines()]
        genotypes_w_missing = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
            [1,  0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0,  1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0,  2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1,  0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1,  1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0,  1, 0, 0, 0, 0, 1, 1, 0, 0, 0]

        ]
        index = 1
        loci = ped_parser.get_loci()

        hetero_freq_tped = [-1, 0.4545, 0.3636, 0.4545, 0.3636, 0.2727, 0.2727]
        for snp in loci:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertAlmostEqual(hetero_freq_tped[index], snp.hetero_freq, places=4)
            index += 1
        self.assertEqual(6, ped_parser.locus_count)
        index = 1
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(genotypes_w_missing[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)



    def testPedBoundaryTPed(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        DataParser.boundary = BoundaryCheck()
        BoundaryCheck.chrom = 2
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        pedigree = [x.split() for x in open(self.tped_filename).readlines()]

        index = 4
        loci = ped_parser.get_loci()
        for snp in loci:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            self.assertAlmostEqual(self.hetero_freq_tped[index], snp.hetero_freq, places=4)

            index += 1
        self.assertEqual(3, ped_parser.locus_count)
        index = 4
        for snp in ped_parser:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            self.assertEqual(pedigree[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))

            index += 1
        self.assertEqual(7, index)

    def testPedSnpBoundaryTPed(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        DataParser.boundary = SnpBoundaryCheck(snps=["rs0001-rs0003"])
        BoundaryCheck.chrom = 1
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        pedigree = [x.split() for x in open(self.tped_filename).readlines()]
        index = 0
        loci = ped_parser.get_loci()
        for snp in loci:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            self.assertAlmostEqual(self.hetero_freq_tped[index], snp.hetero_freq, places=4)
            index += 1
        self.assertEqual(3, ped_parser.locus_count)
        index = 0
        for snp in ped_parser:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            self.assertEqual(pedigree[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))

            index += 1
        self.assertEqual(3, index)

    def testPedSnpBoundary2TPed(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)

        DataParser.boundary = SnpBoundaryCheck(snps=["rs0005-rs0006"])
        BoundaryCheck.chrom = 2
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        pedigree = [x.split() for x in open(self.tped_filename).readlines()]
        index = 4
        loci = ped_parser.get_loci()
        for snp in loci:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            self.assertAlmostEqual(self.hetero_freq_tped[index], snp.hetero_freq, places=4)
            index += 1
        self.assertEqual(2, ped_parser.locus_count)
        index = 4
        for snp in ped_parser:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            self.assertEqual(pedigree[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))

            index += 1
        self.assertEqual(6, index)


    def testPedRegionBoundaryTPed(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)

        DataParser.boundary = SnpBoundaryCheck(snps=["rs0005-rs0006"])
        BoundaryCheck.chrom = 2
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        pedigree = [x.split() for x in open(self.tped_filename).readlines()]
        index = 4
        loci = ped_parser.get_loci()
        for snp in loci:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            index += 1
        index = 4
        for snp in ped_parser:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            self.assertEqual(pedigree[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))

            index += 1
        self.assertEqual(6, index)

    def testPedRegionBoundaryWithExclusionsTPed(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)

        DataParser.boundary = SnpBoundaryCheck(snps=["rs0005-rs0007"])
        DataParser.boundary.LoadExclusions(snps=["rs0007"])
        BoundaryCheck.chrom = 2
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        pedigree = [x.split() for x in open(self.tped_filename).readlines()]
        index = 4
        loci = ped_parser.get_loci()
        for snp in loci:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            index += 1
        index = 4
        for snp in ped_parser:
            self.assertEqual(int(pedigree[index][0]), snp.chr)
            self.assertEqual(int(pedigree[index][3]), snp.pos)
            self.assertEqual(pedigree[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))

            index += 1
        self.assertEqual(6, index)
class TestPedFilesTPedVariedColumns(TestBase):
    def setUp(self):
        super(TestPedFilesTPedVariedColumns, self).setUp()
        self.has_sex = DataParser.has_sex
        self.has_parents = DataParser.has_parents
        self.has_fid = DataParser.has_fid
        self.has_pheno = DataParser.has_pheno
        self.has_liability = DataParser.has_liability

    def tearDown(self):
        super(TestPedFilesTPedVariedColumns, self).tearDown()
        DataParser.has_sex = self.has_sex
        DataParser.has_parents = self.has_parents
        DataParser.has_fid = self.has_fid
        DataParser.has_pheno = self.has_pheno
        DataParser.has_liability = self.has_liability
    # Test to make sure we can load everything
    def testTPedNoFamID(self):
        f = open(self.tfam_filename, "w")
        f.write("""1 0 0 1 0.1
2 0 0 1 0.4
3 0 0 2 1.0
4 0 0 2 0.5
5 0 0 1 0.9
6 0 0 1 1.0
7 0 0 1 0.1
8 0 0 1 0.4
9 0 0 2 1.0
10 0 0 2 0.5
11 0 0 1 0.9
12 0 0 1 1.0""")


        f.close()
        DataParser.has_fid = False
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testTPedNoSex(self):
        f = open(self.tfam_filename, "w")
        f.write("""1 1 0 0 0.1
2 2 0 0 0.4
3 3 0 0 1.0
4 4 0 0 0.5
5 5 0 0 0.9
6 6 0 0 1.0
7 7 0 0 0.1
8 8 0 0 0.4
9 9 0 0 1.0
10 10 0 0 0.5
11 11 0 0 0.9
12 12 0 0 1.0""")


        f.close()
        DataParser.has_sex = False
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testTPedNoParents(self):
        f = open(self.tfam_filename, "w")
        f.write("""1 1 1 0.1
2 2 1 0.4
3 3 2 1.0
4 4 2 0.5
5 5 1 0.9
6 6 1 1.0
7 7 1 0.1
8 8 1 0.4
9 9 2 1.0
10 10 2 0.5
11 11 1 0.9
12 12 1 1.0""")


        f.close()
        DataParser.has_parents = False
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)



    def testTPedNopheno(self):
        f = open(self.tfam_filename, "w")
        f.write("""1 1 0 0 1
2 2 0 0 1
3 3 0 0 2
4 4 0 0 2
5 5 0 0 1
6 6 0 0 1
7 7 0 0 1
8 8 0 0 1
9 9 0 0 2
10 10 0 0 2
11 11 0 0 1
12 12 0 0 1""")
        f.close()
        DataParser.has_pheno = False
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)


    def testTPedLiability(self):
        f = open(self.tfam_filename, "w")
        f.write("""1 1 0 0 1 0.1 1
2 2 0 0 1 0.4 1
3 3 0 0 2 1.0 1
4 4 0 0 2 0.5 1
5 5 0 0 1 0.9 1
6 6 0 0 1 1.0 1
7 7 0 0 1 0.1 1
8 8 0 0 1 0.4 1
9 9 0 0 2 1.0 1
10 10 0 0 2 0.5 1
11 11 0 0 1 0.9 1
12 12 0 0 1 1.0 1""")
        f.close()

        DataParser.has_liability = True
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)




    def testTPedNoFamIDSex(self):
        f = open(self.tfam_filename, "w")
        f.write("""1 0 0 0.1
2 0 0 0.4
3 0 0 1.0
4 0 0 0.5
5 0 0 0.9
6 0 0 1.0
7 0 0 0.1
8 0 0 0.4
9 0 0 1.0
10 0 0 0.5
11 0 0 0.9
12 0 0 1.0""")


        f.close()

        DataParser.has_fid = False
        DataParser.has_sex = False
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testTPedNoParentsPheno(self):
        f = open(self.tfam_filename, "w")
        f.write("""1 1 1
2 2 1
3 3 2
4 4 2
5 5 1
6 6 1
7 7 1
8 8 1
9 9 2
10 10 2
11 11 1
12 12 1""")


        f.close()

        DataParser.has_parents = False
        DataParser.has_pheno   = False
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.tped_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)


class TestTPedFileNegPos(TestBase):
    def testPedNegativePositions(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.misssnp_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.misssnp_tped_filename).readlines()]

        index = 2
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testPedNegativePositionsOtherChrom(self):
        BoundaryCheck.chrom = 2
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.misssnp_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.misssnp_tped_filename).readlines()]

        index = 4
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1

        self.assertEqual(7, index)

    def testPedNegativePositionsLocalChrom(self):
        BoundaryCheck.chrom = 1
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.misssnp_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.misssnp_tped_filename).readlines()]

        index = 2
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(4, index)

    def testPedNegativePosLocalChromMissSNP(self):
        BoundaryCheck.chrom = 1
        DataParser.boundary.LoadExclusions(snps=["rs0004"])

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.misssnp_tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        mapdata = [x.strip().split() for x in open(self.misssnp_tped_filename).readlines()]

        index = 2
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(3, index)

if __name__ == "__main__":
    unittest.main()
