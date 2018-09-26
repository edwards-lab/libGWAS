from data_parser import DataParser
from pheno_covar import PhenoCovar
from boundary import BoundaryCheck
from parsed_locus import ParsedLocus
from exceptions import TooManyAlleles
from exceptions import TooFewAlleles
from bgen_reader import read_bgen
from exceptions import StopIteration
import numpy
import os
from . import ExitIf
from . import BuildReportLine
import sys
import bgen_reader

__copyright__ = "Eric Torstenson"
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

class Parser(DataParser):
    """Parse bgen formatted data

    """
    #: The threshold associated with the .info info column
    info_threshold = 0.4

    def __init__(self, bgen_filename, sample_filename=None):
        """Normal processing will iterate over every locus in each of the files in archive_list

        Client calls should remove anything from archive list that isn't to be parsed"""
        self.bgen_filename = bgen_filename
        self.sample_filename = sample_filename
        self.ind_mask = None
        self.ind_count = -1
        self.bgen_idx = -1
        self.bgen = None            # This is the buffer where we'll store the output from the current file
        self.markers_raw = None     # raw marker details in DataFrame format

        ExitIf("bgen file not found, %s" % (file),
                not os.path.exists(self.bgen_filename))
        if self.sample_filename is not None:
            ExitIf("Sample file, %s, not found. " % (self.sample_filename),
                   not os.path.exists(self.sample_filename))

    def ReportConfiguration(self, file):
        print >> file, BuildReportLine("BGEN FILE", self.bgen_filename)
        if self.sample_file is not None:
            print >> file, BuildReportLine("SAMPLE FILE", self.sample_file)


    def load_family_details(self, pheno_covar):
        """Load family data updating the pheno_covar with  family ids found.

        :param pheno_covar: Phenotype/covariate object
        :return: None
        """

        sample_ids = self.bgen['samples']

        file = open(self.fam_details)
        header = file.readline()
        format = file.readline()        #
        self.file_index = 0

        sample_index = 0
        mask_components = []        # 1s indicate an individual is to be masked out
        for line in file:
            words = line.strip().split()
            indid = PhenoCovar.build_id(words)
            ExitIf(indid != sample_ids[sample_index],
                "Sample id order doesn't match: %s != %s at #%s" % (indid, sample_ids[sample_index], sample_index + 1)
            sample_index += 1
            if DataParser.valid_indid(indid):
                mask_components.append(0)
                pheno_covar.add_subject(indid)
            else:
                mask_components.append(1)
        self.ind_mask = numpy.array(mask_components, dtype=numpy.int8)
        self.ind_count = self.ind_mask.shape[0]
        pheno_covar.freeze_subjects()

    def next_file(self):
        self.file_index += 1
        if self.file_index < len(self.archives):
            self.current_filename = self.archives[self.file_index]
            self.bgen = read_bgen(self.current_filename)

            """self.bgen now holds a dict with the following keys:
            
                    variants => DataFrame with chrom, id, pos, rsid & id
                    genotypes => array (locus x sample x 3xGenoLikelihood)
                    samples => Sample IDs 
                    
                    self.bgen['genotypes'].compute() will return the array in a 
                    format that can be processed
            """
            self.bgen_idx = 0

            if len(self.ids) > 0:
                # We don't care if the IDs match, just that the count is the same
                # If they provide us with mixed up subject IDs, then that's a
                # problem on their end and verifying order of many files could
                # result in a waste of precious computation time
                ExitIf("Inconsistant subject count encountered in file, %s. "% (self.current_filename),
                       len(self.ids) != len(self.bgen['ids']))
        else:
            raise StopIteration

    def open_bgen(self):
        self.bgen = bgen_reader.read_bgen(self.bgen_filename,
                                          sample_file=self.sample_filename,
                                          verbose=False)
        self.markers_raw = bgen['variants']

    def parse_variant(self, index):
        locus = Locus()
        v = self.markers_raw.loc[index]

        locus.chr = v.chrom
        locus.pos = v.pos
        locus.alleles = v.allele_ids.split(",")
        locus.rsid = v.rsid
        if ":" in locus.rsid:
            locus.rsid = locus.rsid.split(":")[0]
        locus.cur_idx = index
        locus.nalleles = locus.nalleles

    # We'll assume that all files that have been associated with this
    # "dataset" are to be considered.
    def load_genotypes(self):
        """Prepares the files for genotype parsing.
        :return: None
        """
        self.open_bgen()

    def get_next_line(self):
        """If we reach the end of the file, we simply open the next, until we \
        run out of archives to process"""
        locus = self.parse_variant(self.bgen_index)
        self.bgen_idx += 1
        return locus, 1.0
        #return line, info, exp_freq

    def populate_iteration(self, iteration):
        """Parse genotypes from the file and iteration with relevant marker \
            details.

        :param iteration: ParseLocus object which is returned per iteration
        :return: True indicates current locus is valid.

        StopIteration is thrown if the marker reaches the end of the file or
        the valid genomic region for analysis.
        """
        global encoding
        snpdata, info = self.get_next_line()

        if info > Parser.info_threshold:
            iteration.chr = snpdata.chr
            iteration.pos = snpdata.pos
            nalleles = snpdata.nalleles
            iteration.rsid = snpdata.rsid
            if DataParser.boundary.TestBoundary(iteration.chr, iteration.pos, iteration.rsid):
                genotypes = numpy.ma.MaskedArray(self.bgen['genotype'][self.bgen_idx - 1].compute(), self.ind_mask).compressed()
                estimate = None
                maf = none
                if encoding == Encoding.Dominant:
                    estimate = genotypes[:, 1] + genotypes[:, 2]
                elif encoding == Encoding.Additive:
                    estimate = genotypes[:, 1] + 2 * genotypes[:, 2]
                    maf = estimate
                elif encoding == Encoding.Recessive:
                    estimate = genotypes[2]
                elif encoding == Encoding.Genotype:
                    estimates = numpy.full_like(genotypes.shape, 2)
                    estimate[genotypes[:, 1] > genotypes[:, 0] and genotypes[:, 1] > genotypes[:, 2]] = 1
                    estimate[genotypes[:, 0] > genotypes[:, 1] and genotypes[:, 0] > genotypes[:, 2]] = 0
                iteration.non_missing_alc = genotypes.shape[0] * 2
                if maf is None:
                    maf = genotypes[:, 1] + 2 * genotypes[:, 2]
                maf = numpy.mean(maf)/2
                iteration.allele_count2 = maf * (iteration.non_missing_alc * 2)
                iteration.effa_freq = maf

                if maf > 0.5:
                    iteration.min_allele_count = iteration.non_missing_alc - iteration.allele_count2
                    iteration.maj_allele_count = iteration.allele_count2
                    maf = 1.0 - maf
                else:
                    iteration.min_allele_count = iteration.allele_count2
                    iteration.maj_allele_count = iteration.non_missing_alc - iteration.allele_count2

                iteration._maf = maf
                iteration.genotype_data = numpy.array(estimate)

                return iteration.maf >= DataParser.min_maf and iteration.maf <= DataParser.max_maf
            else:
                return False
        else:
            return False

    def __iter__(self):
        """Reset the file and begin iteration"""

        return ParsedLocus(self)
