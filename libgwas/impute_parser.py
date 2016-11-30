from . import ExitIf
from . import BuildReportLine
from data_parser import DataParser
from parsed_locus import ParsedLocus
from exceptions import TooManyAlleles
from exceptions import TooFewAlleles
import gzip
import numpy
from exceptions import InvalidSelection
import os

__copyright__ = "Todd Edwards, Chun Li & Eric Torstenson"
__license__ = "GPL3.0"
#     This file is part of libGWAS.
#
#     libGWAS is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     libGWAS is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with MVtest.  If not, see <http://www.gnu.org/licenses/>.

class Encoding(object):
    """Simple enumeration for various model encodings"""
    Additive = 0
    Dominant = 1
    Recessive = 2
    Genotype = 3
    Raw = 4
encoding = Encoding.Additive


def SetEncoding(sval):
    """Sets the encoding variable according to the text passed

    :param sval: text specification for the desired model
    """
    global encoding
    s=sval.lower()
    if s == "additive":
        encoding = Encoding.Additive
    elif s == "dominant":
        encoding = Encoding.Dominant
    elif s == "recessive":
        encoding = Encoding.Recessive
    elif s == "genotype":
        encoding = Encoding.Genotype
    elif s == "raw":
        encoding = Encoding.Raw
    else:
        raise InvalidSelection("Invalid encoding, %s, selected" % (sval))


"""
ISSUES:
* Beyond consideration for MVTest, is it typical to transform these frequencies into genotypes?
* Currently, we will not be filtering on individuals except by explicit removal
* We are assuming that each gzip archive contains all data associated with the loci contained within (i.e. there won't
  be separate files with different subjects inside) (( Todd email sep-6-2014))
* .gen_samples appear to be identical  Any reason to have more than one of them? (even if the user ended up creating a
  bunch in order to toss the job onto a cluster)
* There is no reason to process regions in any order. I'm thinking we'll have a master file and then indices into that
  file and task count to facilitate "parallel" execution
"""

class Parser(DataParser):
    """Parse IMPUTE style output.

        """
    #: the extension associated with the .info files if not using conventions
    info_ext = "info"

    #: The genotype file suffix (of not following convention)
    gen_ext = "gen.gz"

    #: The threshold associated with the .info info column
    info_threshold = 0.4

    def __init__(self, fam_details, archive_list, chroms, info_files=[]):
        """Initialize the structure with the family details file and the list of archives to be parsed

        """
        self.name = None
        ExitIf("Imputed Family file not found, %s" % (fam_details), not os.path.exists(fam_details))
        for file in archive_list:
            ExitIf("Archive file not found, %s" % (file), not os.path.exists(file))
            if self.name is not None:
                self.name = file.split("/")[-1].split(".")[0]
            if len(info_files) == 0:
                info_file = file.replace(Parser.gen_ext, Parser.info_ext)
                ExitIf("Info file not found, %s" % (info_file), not os.path.exists(info_file))
                ExitIf("Info and sample files appear to be same. Is the gen_ext invalid? (%s)" % info_file, info_file == file)
        idx = 0
        for file in info_files:
            ExitIf("Info file not found, %s" % (file), not os.path.exists(file))
            ExitIf("Info filename can't be same as sample filename, %s == %s" % (file, archive_list[idx]), file==archive_list[idx])
            idx += 1

        #: single file containing the subject details (similar to plink's .fam)
        self.fam_details = fam_details

        #: This is only the list of files to be processed
        self.archives = archive_list

        #: array of .info files
        self.info_files = info_files

        #: This will be used to record the opened file used for parsing
        self.current_file = None

        #: This will be used to record the chromosome of the current file
        self.current_chrom = None

        #: This will be used to record the info file associated with quality of SNPs
        self.current_info = None

        #: List of chroms to match files listed in archives
        self.chroms = chroms


    def ReportConfiguration(self, file):
        """
        :param file: Destination for report details
        :return: None
        """
        global encodingpar
        print >> file, BuildReportLine("FAM FILE", self.fam_details)
        print >> file, BuildReportLine("IMPUTE_ARCHIVES", "%s:%s" % (str(self.chroms[0]), self.archives[0]))
        idx = 0
        for arch in self.archives[1:]:
            print >> file, BuildReportLine("", "%s:%s" % (str(self.chroms[idx+1]), arch))
            idx += 1
        print >> file, BuildReportLine("ENCODING", ["Additive", "Dominant", "Recessive", "Genotype", "Raw"][encoding])
        print >> file, BuildReportLine("INFO-EXT", Parser.info_ext)
        print >> file, BuildReportLine("INFO-THRESH", Parser.info_threshold)

    def load_family_details(self, pheno_covar):
        """Load family data updating the pheno_covar with  family ids found.

        :param pheno_covar: Phenotype/covariate object
        :return: None
        """
        file = open(self.fam_details)
        header = file.readline()
        format = file.readline()
        self.file_index = 0

        mask_components = []        # 1s indicate an individual is to be masked out
        for line in file:
            words = line.strip().split()
            indid = ":".join(words[0:2])
            if DataParser.valid_indid(indid):
                mask_components.append(0)
                sex = int(words[5])
                pheno = float(words[6])
                pheno_covar.add_subject(indid, sex, pheno)
            else:
                mask_components.append(1)
        mask_components = numpy.array(mask_components)
        self.ind_mask = numpy.zeros(len(mask_components) * 2, dtype=numpy.int8).reshape(-1, 2)
        self.ind_mask[0:, 0] = mask_components
        self.ind_mask[0:, 1] = mask_components
        self.ind_count = self.ind_mask.shape[0]
        pheno_covar.freeze_subjects()

    def load_genotypes(self):
        """Prepares the files for genotype parsing.

        :return: None
        """


        if self.file_index < len(self.archives):
            self.current_file = self.archives[self.file_index]
            info_filename = self.current_file.replace(Parser.gen_ext, Parser.info_ext)
            if len(self.info_files) > 0:
                info_filename = self.info_files[self.file_index]
            self.info_file = open(info_filename)
            self.info_file.readline()   # Dump the header

            if DataParser.compressed_pedigree:
                self.freq_file = gzip.open("%s" % (self.current_file), 'rb')
            else:
                self.freq_file = open(self.current_file)
            self.current_chrom = self.chroms[self.file_index]
            self.file_index += 1
        else:
            raise StopIteration

    def get_next_line(self):
        """If we reach the end of the file, we simply open the next, until we \
        run out of archives to process"""

        line = self.freq_file.readline().strip().split()
        if len(line) < 1:
            self.load_genotypes()
            line = self.freq_file.readline().strip().split()
        info_line = self.info_file.readline().strip().split()
        info = float(info_line[4])
        exp_freq = float(info_line[3])
        return line, info, exp_freq

    def get_effa_freq(self, genotypes):
        """Returns the effect allele's frequency"""
        return numpy.mean(numpy.array(genotypes))/2

    def populate_iteration(self, iteration):
        """Parse genotypes from the file and iteration with relevant marker \
            details.

        :param iteration: ParseLocus object which is returned per iteration
        :return: True indicates current locus is valid.

        StopIteration is thrown if the marker reaches the end of the file or
        the valid genomic region for analysis.
        """
        global encoding
        line, info, exp_freq = self.get_next_line()

        if info > Parser.info_threshold:
            junk, iteration.rsid, iteration.pos, iteration.major_allele, iteration.minor_allele = line[0:5]
            iteration.chr = self.current_chrom
            iteration.pos = int(iteration.pos)
            if DataParser.boundary.TestBoundary(iteration.chr, iteration.pos, iteration.rsid):
                frequencies = []
                idx = 5
                total_maf = 0.0
                additive = []
                for is_ignored in self.ind_mask[:,0]:
                    if not is_ignored:
                        AA,Aa,aa = [float(x) for x in line[idx:idx+3]]
                        additive.append(Aa+2*aa)
                        if encoding==Encoding.Dominant:
                            estimate = Aa + aa
                        elif encoding==Encoding.Additive:
                            estimate = additive[-1]
                        elif encoding==Encoding.Recessive:
                            estimate = aa
                        elif encoding==Encoding.Genotype:
                            if Aa >= AA and Aa >= aa:
                                estimate = 1
                            elif AA >= Aa and AA >= aa:
                                estimate = 0
                            else:
                                estimate = 2
                        elif encoding==Encoding.Raw:
                            estimate = [AA, Aa, aa]
                        total_maf += numpy.sqrt(aa)
                        frequencies.append(estimate)
                    idx += 3
                iteration.non_missing_alc = len(additive)*2
                maf = numpy.mean(numpy.array(additive))/2
                iteration.allele_count2 = maf * (len(additive) * 2)
                iteration.effa_freq = maf

                if maf > 0.5:
                    iteration.min_allele_count = len(additive)*2 - iteration.allele_count2
                    iteration.maj_allele_count = iteration.allele_count2
                    maf = 1.0 - maf
                else:
                    iteration.min_allele_count = iteration.allele_count2
                    iteration.maj_allele_count = len(additive)*2 - iteration.allele_count2

                iteration._maf = maf
                iteration.genotype_data = numpy.array(frequencies)

                return iteration.maf >= DataParser.min_maf and iteration.maf <= DataParser.max_maf
            else:
                return False
        else:
            return False

    def __iter__(self):
        """Reset the file and begin iteration"""

        return ParsedLocus(self)
