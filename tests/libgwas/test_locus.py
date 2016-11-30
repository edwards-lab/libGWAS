import test_transped_parser
from libgwas.data_parser import DataParser
from libgwas.pheno_covar import PhenoCovar
from libgwas.transposed_pedigree_parser import Parser as TransposedPedigreeParser
from libgwas.locus import Locus


class TestLocusBasics(test_transped_parser.TestBase):
    # Test to make sure we can load everything
    def testAllelesIteration(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()


        index = 0
        for snp in ped_parser:
            self.assertEqual(self.tped1_alleles[index][1], snp.minor_allele)
            self.assertEqual(self.tped1_alleles[index][0], snp.major_allele)

            index += 1
        self.assertEqual(7, index)

    # There was an error where the alleles array is updated despite the attempt at storing actual copies of the
    # Loci in the locus vector. This reflects that error.
    def testAllelesInLoci(self):
        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()


        index = 0
        for snp in ped_parser.get_loci():
            self.assertEqual(self.tped1_alleles[index][1], snp.minor_allele)
            self.assertEqual(self.tped1_alleles[index][0], snp.major_allele)

            index += 1
        self.assertEqual(7, index)

    def testLocusRelation(self):
        l1 = Locus()
        l1.chr = 1
        l1.pos = 100

        l2 = Locus()
        l2.chr = 1
        l2.pos = 101

        l3 = Locus()
        l3.chr = 2
        l3.pos = 50

        l4 = Locus()
        l4.chr = 1
        l4.pos = 100

        self.assertTrue(l1 == l1)
        self.assertTrue(l1 == l4)
        self.assertTrue(l1 < l2)
        self.assertTrue(l2 > l1)
        self.assertTrue(l1 < l3)
        self.assertTrue(l3 > l1)
        self.assertFalse(l1 == l2)
        self.assertFalse(l2 < l1)

        s = set([l1, l2, l3])
        self.assertTrue(l4 in s)

        d = {l1:"One", l2:"Two", l3:"Three"}
        self.assertTrue(l4 in d)
        self.assertTrue(d[l4] == "One")

