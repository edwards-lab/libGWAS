#!/usr/bin/env python
import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../../")
    sys.path.insert(0, "../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")
import os
import unittest
import numpy
from libgwas.pheno_covar import PhenoCovar
from libgwas.exceptions import MalformedInputFile
from libgwas.exceptions import InvalidSelection
from libgwas.exceptions import NoMatchedPhenoCovars
from libgwas.data_parser import DataParser
import libgwas.standardizer

class TestBase(unittest.TestCase):
    def setUp(self):
        self.pcsac = PhenoCovar.sex_as_covariate
        self.filenames = self.WriteTestFiles()

        self.no_header = open(self.filenames[0])

        self.header = open(self.filenames[1])
        self.ped = [l.strip() for l in open(self.filenames[2]).readlines()]
        self.nh_filename = self.filenames[0]        # No header
        self.h_filename = self.filenames[1]         # Header
        self.ped_filename = self.filenames[2]       # Pedigree
        self.mch_filename = self.filenames[4]       # Multicolumn with header
        self.mc_filename = self.filenames[3]        # multicolumn without header

        self.phenotypes = [[0.1, 1.0, 0.5],
             [0.2, 0.5, 1.0],
             [0.3, 0.6, 0.1],
             [0.4, 0.5, 0.5],
             [0.5, 1.0, 1.0],
             [0.6, 0.1, 0.2]]
        self.sex = [1,1,2,2,1,1]
        self.has_pheno = DataParser.has_pheno
        DataParser.has_pheno = False
        self.standardizer = libgwas.standardizer.get_standardizer()
        libgwas.standardizer.set_standardizer(libgwas.standardizer.NoStandardization)




    def tearDown(self):
        for file in self.filenames:
            os.remove(file)
        PhenoCovar.sex_as_covariate = self.pcsac
        DataParser.has_pheno = self.has_pheno
        libgwas.standardizer.set_standardizer(self.standardizer)

    def WriteTestFiles(self, prefix = "__test_pheno"):
        filenames = []
        filename = "%s_sc.txt" % (prefix)
        f = open(filename, "w")
        f.write("""Fam1\tInd1\t0.9
Fam2\tInd2\t1.0
Fam3\tInd3\t0.4
Fam4\tInd4\t0.8
Fam4\tInd5\t1
Fam4\tInd6\t0.1""")
        f.close()
        filenames.append(filename)

        filename = "%s_sch.txt" % (prefix)
        f = open(filename, "w")
        f.write("""fid\tIID\tBMI
Fam1\tInd1\t0.9
Fam2\tInd2\t1.0
Fam3\tInd3\t0.4
Fam4\tInd4\t0.8
Fam4\tInd5\t1
Fam4\tInd6\t0.1""")
        f.close()
        filenames.append(filename)

        filename = "%s.ped" % (prefix)
        f = open(filename, "w")
        f.write("""Fam1\tInd1\t1\t0.1
Fam2\tInd2\t1\t0.4
Fam3\tInd3\t2\t1.0
Fam4\tInd4\t2\t0.5
Fam4\tInd5\t1\t0.9
Fam4\tInd6\t1\t1.0""")
        f.close()
        filenames.append(filename)



        filename = "%s_mc.txt" % (prefix)
        f = open(filename, "w")
        f.write("""Fam1\tInd1\t0.1\t1.0\t0.5
Fam2\tInd2\t0.2\t0.5\t1.0
Fam3\tInd3\t0.3\t0.6\t0.1
Fam4\tInd4\t0.4\t0.5\t0.5
Fam4\tInd5\t0.5\t1.0\t1.0
Fam4\tInd6\t0.6\t0.1\t0.2""")
        f.close()
        filenames.append(filename)



        filename = "%s_mch.txt" % (prefix)
        f = open(filename, "w")
        f.write("""FID\tIID\tBMI\tIBM\tMSA
Fam1\tInd1\t0.1\t1.0\t0.5
Fam2\tInd2\t0.2\t0.5\t1.0
Fam3\tInd3\t0.3\t0.6\t0.1
Fam4\tInd4\t0.4\t0.5\t0.5
Fam4\tInd5\t0.5\t1.0\t1.0
Fam4\tInd6\t0.6\t0.1\t0.2""")
        f.close()
        filenames.append(filename)



        # For verification that we handle mismatched IDs
        filename = "%s-mm.ped" % (prefix)
        f = open(filename, "w")
        f.write("""Fam1\tInd1\t1\t0.1
Fam2\tInd2\t1\t0.4
Fam3\tInd3\t2\t1.0
Fam4\tInd4\t2\t0.5
Fam4\tInd5\t1\t0.9
Fam4\tInd6\t1\t1.0
Fam4\tInd7\t2\t-9
Fam9\tInd1\t1\t0.9""")
        f.close()
        filenames.append(filename)



        filename = "%s_mm.txt" % (prefix)
        f = open(filename, "w")
        f.write("""FID\tIID\tBMI
Fam0\tInd2\t-9
Fam1\tInd1\t0.9
Fam2\tInd2\t1.0
Fam3\tInd3\t0.4
Fam4\tInd4\t0.8
Fam4\tInd5\t1
Fam4\tInd6\t0.1""")
        f.close()
        filenames.append(filename)


        filename = "%s_miss.txt" % (prefix)
        f = open(filename, "w")
        f.write("""FID\tIID\tBMI\tIBM\tMSA
Fam1\tInd1\t-9\t1.0\t0.5
Fam2\tInd2\t0.2\t-9\t1.0
Fam3\tInd3\t0.3\t0.6\t-9
Fam4\tInd4\t0.4\t0.5\t0.5
Fam4\tInd5\t0.5\t1.0\t1.0
Fam4\tInd6\t0.6\t0.1\t0.2""")
        f.close()
        filenames.append(filename)

        filename = "%s_mch.sample" % (prefix)
        f = open(filename, "w")
        f.write("""FID\tIID\tBMI\tIBM\tMSA
0\t0\tC\tC\tC
Fam1\tInd1\t0.1\t1.0\t0.5
Fam2\tInd2\t0.2\t0.5\t1.0
Fam3\tInd3\t0.3\t0.6\t0.1
Fam4\tInd4\t0.4\t0.5\t0.5
Fam4\tInd5\t0.5\t1.0\t1.0
Fam4\tInd6\t0.6\t0.1\t0.2""")
        f.close()
        filenames.append(filename)

        filename = "%s_sch.sample" % (prefix)
        f = open(filename, "w")
        f.write("""fid\tIID\tBMI
0\t0\tC
Fam1\tInd1\t0.9
Fam2\tInd2\t1.0
Fam3\tInd3\t0.4
Fam4\tInd4\t0.8
Fam4\tInd5\t1
Fam4\tInd6\t0.1""")
        f.close()
        filenames.append(filename)

        return filenames


def load_pedigree(pc, ped):
    for line in ped:
        fam, ind, sex, ph = line.split()
        pc.add_subject("%s:%s" % (fam, ind), sex=int(sex), phenotype=float(ph))
    pc.freeze_subjects()
class TestIteration(TestBase):
    def setUp(self):
        super(TestIteration, self).setUp()
        # This class deals with a different set of files for header/pedigree information
        # We want some individuals with data missing on both sides
        self.header = open(self.filenames[4])
        self.ped = [l.strip() for l in open(self.filenames[5]).readlines()]

        self.phenotypes = [[ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                           [ 1.0, 0.5, 0.6, 0.5, 1.0, 0.1],
                           [ 0.5, 1.0, 0.1, 0.5, 1.0, 0.2]]

    def testEmptyIterator(self):
        pc = PhenoCovar()
        count = 0
        for test in pc:
            count += 1

        self.assertEqual(0, count)


    def testBasicIteration(self):
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        orig_pheno = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, 0.9]
        sex   = [1, 1, 2, 2, 1, 1, 1]

        count = 0

        for test in pc:
            pheno, covars, nonmissing = test.get_variables()
            self.assertEqual(8, len(nonmissing))
            self.assertEqual(7, numpy.sum(nonmissing))
            self.assertEqual(1, len(covars))

            self.assertEqual(sex, list(covars[count]))
            self.assertEqual("SEX", test.get_covariate_name(count))
            self.assertEqual(7, len(pheno))
            self.assertEqual("Pheno-1", test.get_phenotype_name())
            for i in range(0, len(pheno)):
                self.assertAlmostEqual(orig_pheno[i], pheno[i])
            count += 1
        self.assertEqual(1, count)


    def testBasicWithMask(self):
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.individual_mask = [0, 0, 1, 0, 0, 0, 1, 0]

        orig_pheno = [0.1, 0.4, 0.5, 0.9, 1.0, 0.9]
        sex   = [1, 1, 2, 1, 1, 1]

        count = 0
        for test in pc:
            pheno, covars, nonmissing = test.get_variables(numpy.array(pc.individual_mask, dtype=bool))
            self.assertEqual(6, numpy.sum(nonmissing))
            self.assertEqual(1, len(covars))

            self.assertEqual(sex, list(covars[count]))
            self.assertEqual("SEX", test.get_covariate_name(count))
            self.assertEqual(6, len(pheno))
            self.assertEqual("Pheno-1", test.get_phenotype_name())
            for idx in range(0, len(pheno)):
                self.assertEqual(orig_pheno[idx], pheno[idx])
            count += 1
        self.assertEqual(1, count)


    def testBasicWithCovar(self):
        PhenoCovar.sex_as_covariate = False
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.individual_mask = [0, 0, 0, 0, 0, 0, 1, 1]
        pc.load_covarfile(self.header, names=["BMI", "MSA"])
        pc.freeze_subjects()

        pheno = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0]
        sex   = [1, 1, 2, 1, 1, 1]

        # First, test without sex as covariate
        count = 0
        for test in pc:
            test_pheno, covars, nonmissing = test.get_variables()
            self.assertEqual(6, numpy.sum(nonmissing))
            self.assertEqual(2, len(covars))
            for phenotype in self.phenotypes[0:2]:
                self.assertEqual(self.phenotypes[0], list(covars[0]))
                self.assertEqual(self.phenotypes[2], list(covars[1]))
                self.assertEqual("BMI", test.get_covariate_name(0))
                self.assertEqual("MSA", test.get_covariate_name(1))
                self.assertEqual(6, len(pheno))
                self.assertEqual("Pheno-1", test.get_phenotype_name())

                for i in range(0, len(pheno)):
                    self.assertAlmostEqual(pheno[i], test_pheno[i])
            count += 1
        self.assertEqual(1, count)


        PhenoCovar.sex_as_covariate = True


class TestMismatchedIDs(TestBase):
    def setUp(self):
        super(TestMismatchedIDs, self).setUp()
        # This class deals with a different set of files for header/pedigree information
        self.header = open(self.filenames[6])
        self.ped = [l.strip() for l in open(self.filenames[5]).readlines()]
        self.pheno = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, PhenoCovar.missing_encoding, 0.9]
        self.sex   = [1, 1, 2, 2, 1, 1, 2, 1]

    def testPhenoHeader(self):
        """Phenotype with header and covariate.

        We are just making sure everything looks reasonable when no complications are involved"""

        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(8, len(pc.phenotype_data[0]))
        self.assertEqual(1, len(pc.phenotype_names))
        self.assertEqual("Pheno-1", pc.phenotype_names[0])

        for idx in xrange(0, len(self.pheno)):
            self.assertAlmostEqual(self.pheno[idx], pc.phenotype_data[0][idx])
            self.assertEqual(self.sex[idx], pc.covariate_data[0][idx])

        pc.load_phenofile(self.header)
        self.assertEqual("BMI", pc.phenotype_names[0])
        pheno =  [0.9, 1.0, 0.4, 0.8, 1, 0.1, PhenoCovar.missing_encoding, PhenoCovar.missing_encoding]
        for idx in xrange(0, len(pheno)):
            self.assertAlmostEqual(pheno[idx], pc.phenotype_data[0][idx])


    def testCovarHeader(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(8, len(pc.phenotype_data[0]))
        self.assertEqual(1, len(pc.phenotype_names))
        self.assertEqual("Pheno-1", pc.phenotype_names[0])

        for idx in xrange(0, len(self.pheno)):
            self.assertAlmostEqual(self.pheno[idx], pc.phenotype_data[0][idx])
            self.assertEqual(self.sex[idx], pc.covariate_data[0][idx])

        pc.load_covarfile(self.header)
        self.assertEqual("BMI", pc.covariate_labels[1])
        self.assertEqual("SEX", pc.covariate_labels[0])
        covar =  [0.9, 1.0, 0.4, 0.8, 1, 0.1, PhenoCovar.missing_encoding, PhenoCovar.missing_encoding]
        for idx in xrange(0, len(covar)):
            self.assertAlmostEqual(self.sex[idx], pc.covariate_data[0][idx])
            self.assertAlmostEqual(covar[idx], pc.covariate_data[1][idx])
            self.assertAlmostEqual(self.pheno[idx], pc.phenotype_data[0][idx])

    def testCovarMissingAll(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = True

        prefix = "__test_pheno"
        filename = "%s_miss.txt" % (prefix)
        f = open(filename, "w")
        f.write("""FID\tIID\tBMI\tIBM\tMSA
F1\tI1\t-9\t1.0\t0.5
F2\tI2\t0.2\t-9\t1.0
F3\tI3\t0.3\t0.6\t-9
F4\tI4\t0.4\t0.5\t0.5
F4\tI5\t0.5\t1.0\t1.0
F4\tI6\t0.6\t0.1\t0.2""")
        f.close()
        file = open(filename)

        with self.assertRaises(NoMatchedPhenoCovars):
            pc = PhenoCovar()
            load_pedigree(pc, self.ped)
            pc.load_covarfile(file, indices=[1])



class TestPedigreePopulating(TestBase):
    def test_basic_population(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        #sex = [1, 1, 2, 2, 1, 1]

        for line in self.ped:
            fam, ind, sex, ph = line.split()
            pc.add_subject("%s:%s" % (fam, ind), sex=int(sex), phenotype=float(ph))

        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        i = 0

        for line in self.ped:
            fam, ind, sex, ph = line.split()
            iid = "%s:%s" % (fam, ind)
            pdata = pc.pedigree_data[iid]
            self.assertEqual(pdata, i)
            self.assertAlmostEqual(float(ph), pc.phenotype_data[0][i])
            self.assertEqual(int(sex), pc.covariate_data[0][i])

            i+=1


        # Indicate that we do not want to use sex as a covariate
        PhenoCovar.sex_as_covariate = False
        newpc = PhenoCovar()

        for line in self.ped:
            fam, ind, sex, ph = line.split()
            newpc.add_subject("%s:%s" % (fam, ind), sex=int(sex), phenotype=float(ph))

        # Test that sex wasn't loaded as a covariate due to the setting of PhenoCovar.sex_as_covariate
        self.assertEqual(0, len(newpc.covariate_data))
        self.assertEqual(1, len(newpc.phenotype_data))
        self.assertEqual(6, len(newpc.pedigree_data))
        self.assertEqual(6, len(newpc.phenotype_data[0]))

        i = 0

        for line in self.ped:
            fam, ind, sex, ph = line.split()
            iid = "%s:%s" % (fam, ind)
            pdata = pc.pedigree_data[iid]
            self.assertEqual(pdata, i)
            self.assertAlmostEqual(float(ph), pc.phenotype_data[0][i])

            i+=1


    def test_pheno_with_header(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_phenofile(open(self.filenames[4]), indices=[1,2])
        self.assertEqual(2, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(6, len(pc.phenotype_data[1]))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(["BMI", "IBM"], pc.phenotype_names)



        index = 0
        for p1, p2, p3 in self.phenotypes:
            self.assertEqual(p1, pc.phenotype_data[0][index])
            self.assertEqual(p2, pc.phenotype_data[1][index])
            index += 1
        pc.load_phenofile(open(self.filenames[4]), indices=[1,2,3])
        self.assertEqual(3, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(6, len(pc.phenotype_data[1]))
        self.assertEqual(6, len(pc.phenotype_data[2]))
        self.assertEqual(["BMI", "IBM", "MSA"], pc.phenotype_names)

        index = 0
        for p1, p2, p3 in self.phenotypes:
            self.assertEqual(p1, pc.phenotype_data[0][index])
            self.assertEqual(p2, pc.phenotype_data[1][index])
            self.assertEqual(p3, pc.phenotype_data[2][index])
            index += 1

    def test_sample_pheno_with_header(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_phenofile(open(self.filenames[8]), indices=[1,2], sample_file=True)
        self.assertEqual(2, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(6, len(pc.phenotype_data[1]))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(["BMI", "IBM"], pc.phenotype_names)



        index = 0
        for p1, p2, p3 in self.phenotypes:
            self.assertEqual(p1, pc.phenotype_data[0][index])
            self.assertEqual(p2, pc.phenotype_data[1][index])
            index += 1
        pc.load_phenofile(open(self.filenames[8]), indices=[1,2,3], sample_file=True)
        self.assertEqual(3, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(6, len(pc.phenotype_data[1]))
        self.assertEqual(6, len(pc.phenotype_data[2]))
        self.assertEqual(["BMI", "IBM", "MSA"], pc.phenotype_names)

        index = 0
        for p1, p2, p3 in self.phenotypes:
            self.assertEqual(p1, pc.phenotype_data[0][index])
            self.assertEqual(p2, pc.phenotype_data[1][index])
            self.assertEqual(p3, pc.phenotype_data[2][index])
            index += 1


class TestPhenoWithMissing(TestBase):
    def test_pheno_with_missing(self):
        PhenoCovar.sex_as_covariate = False

        phenotypes = [[-9, 1.0, 0.5],
             [0.2, -9, 1.0],
             [0.3, 0.6, -9],
             [0.4, 0.5, 0.5],
             [0.5, 1.0, 1.0],
             [0.6, 0.1, 0.2]]

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_phenofile(open(self.filenames[7]), indices=[1])
        pc.load_covarfile(open(self.filenames[7]), indices=[2,3])

        missing = [True,  True,  True, False, False, False]
        index = 0
        for y in pc:
            pheno, covars, nonmissing = y.get_variables()
            self.assertEqual(3, numpy.sum(nonmissing))
            self.assertEqual(missing, list(numpy.invert(nonmissing)))
            #non_missing = numpy.ones((6), dtype=bool)
            #self.assertEqual(missing, list(y.missing))
            #(pheno, covariates) = y.get_variables(non_missing)




class TestCovarSingleCol(TestBase):

    def test_exceed_column_count(self):
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        self.assertRaises(InvalidSelection, pc.load_covarfile, self.header, indices=[12])
        self.assertRaises(InvalidSelection, pc.load_covarfile, self.no_header, indices=[12])

    def test_covar_no_header(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = False
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_covarfile(self.no_header)
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(1, len(pc.phenotype_names))
        self.assertEqual("Cov-1", pc.covariate_labels[0])


        covar_values = [0.9, 1.0, 0.4, 0.8, 1, 0.1]
        for idx in xrange(0, len(covar_values)):
            self.assertAlmostEqual(covar_values[idx], pc.covariate_data[0][idx])

        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_covarfile(self.no_header)
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(2, len(pc.covariate_data))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(6, len(pc.covariate_data[1]))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(2, len(pc.covariate_labels))
        self.assertEqual("SEX", pc.covariate_labels[0])
        self.assertEqual("Cov-1", pc.covariate_labels[1])

        sex = [1,1,2,2,1,1]
        for idx in xrange(0, len(covar_values)):
            self.assertAlmostEqual(sex[idx], pc.covariate_data[0][idx])
            self.assertAlmostEqual(covar_values[idx], pc.covariate_data[1][idx])

    def test_covar_header(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = False
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        pc.load_covarfile(self.header)

        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(1, len(pc.covariate_labels))
        self.assertEqual("BMI", pc.covariate_labels[0])

        covariate_values = [0.9, 1.0, 0.4, 0.8, 1, 0.1]
        for idx in xrange(0, len(covariate_values)):
            self.assertAlmostEqual(covariate_values[idx], pc.covariate_data[0][idx])


        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        pc.load_covarfile(self.header)

        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(2, len(pc.covariate_labels))
        self.assertEqual("SEX", pc.covariate_labels[0])
        self.assertEqual("BMI", pc.covariate_labels[1])

        sex = [1,1,2,2,1,1]
        covariate_values = [0.9, 1.0, 0.4, 0.8, 1, 0.1]

        for idx in xrange(0, len(covariate_values)):
            self.assertEqual(sex[idx], pc.covariate_data[0][idx])
            self.assertAlmostEqual(covariate_values[idx], pc.covariate_data[1][idx])

    def test_sample_covar_header(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = False
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        file = open(self.filenames[9])

        pc.load_covarfile(file, sample_file=True)

        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(1, len(pc.covariate_labels))
        self.assertEqual("BMI", pc.covariate_labels[0])

        covariate_values = [0.9, 1.0, 0.4, 0.8, 1, 0.1]
        for idx in xrange(0, len(covariate_values)):
            self.assertAlmostEqual(covariate_values[idx], pc.covariate_data[0][idx])


        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        pc.load_covarfile(file, sample_file=True)

        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(2, len(pc.covariate_labels))
        self.assertEqual("SEX", pc.covariate_labels[0])
        self.assertEqual("BMI", pc.covariate_labels[1])

        sex = [1,1,2,2,1,1]
        covariate_values = [0.9, 1.0, 0.4, 0.8, 1, 0.1]

        for idx in xrange(0, len(covariate_values)):
            self.assertEqual(sex[idx], pc.covariate_data[0][idx])
            self.assertAlmostEqual(covariate_values[idx], pc.covariate_data[1][idx])

    def test_covar_header_with_empty_names(self):
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        # Application was feeding in empty string for names and excepting
        pc.load_covarfile(self.header, names=' ')
        self.assertEqual("BMI", pc.covariate_labels[0])


    def test_covar_header_with_string_indices(self):
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        # Application was feeding in non-integers  for indices and excepting
        pc.load_covarfile(self.header, indices='1')
        self.assertEqual("BMI", pc.covariate_labels[0])

class TestPhenoSingleCol(TestBase):
    def test_exceed_column_count(self):
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        self.assertRaises(InvalidSelection, pc.load_phenofile, self.header, indices=[12])
        self.assertRaises(InvalidSelection, pc.load_phenofile, self.no_header, indices=[12])

    def test_pheno_no_header(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = False
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        pc.load_phenofile(self.no_header)
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(0, len(pc.covariate_data))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(1, len(pc.phenotype_names))
        self.assertEqual("Pheno-1", pc.phenotype_names[0])

        phenotype_values = [0.9, 1.0, 0.4, 0.8, 1, 0.1]
        for idx in xrange(0, len(phenotype_values)):
            self.assertAlmostEqual(phenotype_values[idx], pc.phenotype_data[0][idx])

    def test_pheno_header(self):
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        pc.load_phenofile(self.header)

        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(1, len(pc.phenotype_names))
        self.assertEqual("BMI", pc.phenotype_names[0])

        phenotype_values = [0.9, 1.0, 0.4, 0.8, 1, 0.1]
        for idx in xrange(0, len(phenotype_values)):
            self.assertAlmostEqual(phenotype_values[idx], pc.phenotype_data[0][idx])

    def test_pheno_header_with_empty_names(self):
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        # Application was feeding in empty string for names and excepting
        pc.load_phenofile(self.header, names=' ')


    def test_pheno_header_with_string_indices(self):
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        # Application was feeding in empty string for names and excepting
        pc.load_phenofile(self.header, indices='1')



    def test_invalid_input(self):
        pc = PhenoCovar()

        with self.assertRaises(InvalidSelection):
            pc.load_phenofile(self.header, indices=[2])
        with self.assertRaises(InvalidSelection):
            pc.load_phenofile(self.header, indices=[1, 2, 3])

        with self.assertRaises(MalformedInputFile):
            pc.load_phenofile(self.no_header, names=["BMI"])
        with self.assertRaises(MalformedInputFile):
            pc.load_phenofile(self.no_header, names=["BMI", "AGE"])

        # Make sure that we throw an exception when we provide a name that isn't in the header
        with self.assertRaises(InvalidSelection):
            pc.load_phenofile(self.header, names=["AGE"])


class TestCovarMultiCol(TestBase):
    def testCovarNoHeader(self):
        self.ped = [l.strip() for l in open(self.filenames[2]).readlines()]
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_covarfile(open(self.filenames[3]), indices=[2,3])
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(3, len(pc.covariate_data))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(6, len(pc.covariate_data[1]))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(["Pheno-1"], pc.phenotype_names)
        self.assertEqual(["SEX", "Cov-2", "Cov-3"], pc.covariate_labels)

        sex = [1,1,2,2,1,1]
        covariates = [[0.1, 1.0, 0.5],
             [0.2, 0.5, 1.0],
             [0.3, 0.6, 0.1],
             [0.4, 0.5, 0.5],
             [0.5, 1.0, 1.0],
             [0.6, 0.1, 0.2]]

        index = 0
        for p1, p2, p3 in covariates:
            self.assertEqual(sex[index], pc.covariate_data[0][index])
            self.assertEqual(p2, pc.covariate_data[1][index])
            self.assertEqual(p3, pc.covariate_data[2][index])
            index += 1

        PhenoCovar.sex_as_covariate = False
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)

        pc.load_covarfile(open(self.filenames[3]), indices=[1,2,3])
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(3, len(pc.covariate_data))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(6, len(pc.covariate_data[1]))
        self.assertEqual(6, len(pc.covariate_data[2]))
        self.assertEqual(["Cov-1", "Cov-2", "Cov-3"], pc.covariate_labels)

        index = 0
        for p1, p2, p3 in covariates:
            self.assertEqual(p1, pc.covariate_data[0][index])
            self.assertEqual(p2, pc.covariate_data[1][index])
            self.assertEqual(p3, pc.covariate_data[2][index])
            index += 1
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_covarfile(open(self.filenames[3]), indices=[2,3])
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(3, len(pc.covariate_data))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(6, len(pc.covariate_data[1]))
        self.assertEqual(6, len(pc.covariate_data[2]))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(["Pheno-1"], pc.phenotype_names)
        self.assertEqual(["SEX", "Cov-2", "Cov-3"], pc.covariate_labels)

        index = 0
        for p1, p2, p3 in covariates:
            self.assertEqual(sex[index], pc.covariate_data[0][index])
            self.assertEqual(p2, pc.covariate_data[1][index])
            self.assertEqual(p3, pc.covariate_data[2][index])
            index += 1

    def testCovarNoSex(self):
        self.ped = [l.strip() for l in open(self.filenames[2]).readlines()]

        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = False

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_covarfile(open(self.filenames[4]), indices=[2,3])
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(2, len(pc.covariate_data))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(6, len(pc.covariate_data[1]))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(["Pheno-1"], pc.phenotype_names)
        self.assertEqual(["IBM", "MSA"], pc.covariate_labels)


        index = 0
        for p1, p2, p3 in self.phenotypes:
            self.assertEqual(p2, pc.covariate_data[0][index])
            self.assertEqual(p3, pc.covariate_data[1][index])
            index += 1

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_covarfile(open(self.filenames[4]), indices=[1,2,3])
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(3, len(pc.covariate_data))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(6, len(pc.covariate_data[1]))
        self.assertEqual(6, len(pc.covariate_data[2]))
        self.assertEqual(["BMI", "IBM", "MSA"], pc.covariate_labels)

        index = 0
        for p1, p2, p3 in self.phenotypes:
            self.assertEqual(p1, pc.covariate_data[0][index])
            self.assertEqual(p2, pc.covariate_data[1][index])
            self.assertEqual(p3, pc.covariate_data[2][index])
            index += 1
        PhenoCovar.sex_as_covariate = True
        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_covarfile(open(self.filenames[4]), indices=[2,3])
        self.assertEqual(1, len(pc.phenotype_data))
        self.assertEqual(3, len(pc.covariate_data))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(6, len(pc.covariate_data[1]))
        self.assertEqual(6, len(pc.covariate_data[2]))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(["Pheno-1"], pc.phenotype_names)
        self.assertEqual(["SEX", "IBM", "MSA"], pc.covariate_labels)

        index = 0
        for p1, p2, p3 in self.phenotypes:
            self.assertEqual(self.sex[index], pc.covariate_data[0][index])
            self.assertEqual(p2, pc.covariate_data[1][index])
            self.assertEqual(p3, pc.covariate_data[2][index])
            index += 1


class TestPhenoMultiCol(TestBase):
    def testMultiPhenoHeader(self):
        self.ped = [l.strip() for l in open(self.filenames[2]).readlines()]
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_phenofile(open(self.filenames[4]), indices=[2,3])
        self.assertEqual(2, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(6, len(pc.phenotype_data[1]))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(["IBM", "MSA"], pc.phenotype_names)
        self.assertEqual(["SEX"], pc.covariate_labels)

        sex = [1,1,2,2,1,1]
        covariates = [[0.1, 1.0, 0.5],
             [0.2, 0.5, 1.0],
             [0.3, 0.6, 0.1],
             [0.4, 0.5, 0.5],
             [0.5, 1.0, 1.0],
             [0.6, 0.1, 0.2]]

        index = 0
        for p1, p2, p3 in covariates:
            self.assertEqual(sex[index], pc.covariate_data[0][index])
            self.assertEqual(p2, pc.phenotype_data[0][index])
            self.assertEqual(p3, pc.phenotype_data[1][index])
            index += 1
    def testPhenoNoHeader(self):
        # Indicate that we want to use sex as a covariate
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        load_pedigree(pc, self.ped)
        pc.load_phenofile(open(self.filenames[3]), indices=[2,3])
        self.assertEqual(2, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(6, len(pc.phenotype_data[1]))
        self.assertEqual(6, len(pc.covariate_data[0]))
        self.assertEqual(["Pheno-2", "Pheno-3"], pc.phenotype_names)

        index = 0
        for p1, p2, p3 in self.phenotypes:
            self.assertEqual(p2, pc.phenotype_data[0][index])
            self.assertEqual(p3, pc.phenotype_data[1][index])
            index += 1
        pc.load_phenofile(open(self.filenames[3]), indices=[1,2,3])
        self.assertEqual(3, len(pc.phenotype_data))
        self.assertEqual(1, len(pc.covariate_data))
        self.assertEqual(6, len(pc.phenotype_data[0]))
        self.assertEqual(6, len(pc.phenotype_data[1]))
        self.assertEqual(6, len(pc.phenotype_data[2]))
        self.assertEqual(["Pheno-1", "Pheno-2", "Pheno-3"], pc.phenotype_names)

        index = 0
        for p1, p2, p3 in self.phenotypes:
            self.assertEqual(p1, pc.phenotype_data[0][index])
            self.assertEqual(p2, pc.phenotype_data[1][index])
            self.assertEqual(p3, pc.phenotype_data[2][index])
            index += 1


if __name__ == "__main__":
    unittest.main()
