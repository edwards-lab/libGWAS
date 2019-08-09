from .data_parser import DataParser
from .pheno_covar import PhenoCovar
from .boundary import BoundaryCheck
from .parsed_locus import ParsedLocus
from .exceptions import TooManyAlleles
from .exceptions import TooFewAlleles
from bgen_reader import read_bgen
from .impute_parser import gen_dosage_extraction
# from exceptions import StopIteration
import numpy
import os
from . import ExitIf
from . import BuildReportLine
#from . import timer
import sys
import bgen_reader
from . import impute_parser
from .locus import Locus
import libgwas
import logging

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

encoding = impute_parser.Encoding.Additive

class Parser(DataParser):
    """Parse bgen formatted data

    ** Note ** qctool and bgen_reader both assume single field IDs. qctool will ignore column 2 when converting
    gen files, so we'll assume that's the standard format and corresponding applications should default to
    PhenoCovar.PhenoIdFormat.FID, otherwise, it's likely that the ids from the phenotype parser will not
    agree with the what's recognized by the bgen file. If there are no sample IDs embedded within the

    """
    #: The threshold associated with the .info info column
    info_threshold = 0.4
    het_threshold = 0.8          # Jackie's suggestion. We can switch back to 1.0.

    # It is possible that there are no chromosomes in the locus details (such as gen files converted to bgen). In
    # these cases, we assume that the user will be working on a file that has only one chromosome and that
    # chromosome was specified by the application
    default_chromosome = -1

    def __init__(self, bgen_filename, sample_filename=None, meta_filename=None):
        """Support is present only for a single .bgen file (and possibly corresponding sample file)

        If sample file is not present, bgen file should have sample IDs baked into
        (which is part of the format).

        Currently, support for the metadata doesn't exist, but files are present
        as placeholders. """
        self.bgen_filename = bgen_filename
        self.sample_filename = sample_filename
        self.meta_filename = meta_filename
        self.ind_mask = None
        self.geno_mask = None       # this will be 3x in order to permit masking of rows within the 3xN vector of gt probs
        self.ind_count = -1
        self.bgen_idx = -1
        self.bgen_start_idx = 0     # We'll mark this for the first locus to be analyzed

        self.bgen = None            # This is the buffer where we'll store the output from the current file
        self.markers_raw = None     # raw marker details in DataFrame format
        # We don't want to bother checking for this until we know which subjects
        # to exclude
        self.max_missing_geno = None
        self.sample_ids = None
        self.alt_not_missing = None

        ExitIf("bgen file not found, %s" % (self.bgen_filename),
                not os.path.exists(self.bgen_filename))
        if self.sample_filename is not None:
            ExitIf("Sample file, %s, not found. " % (self.sample_filename),
                   not os.path.exists(self.sample_filename))
        self.open_bgen()


    def ReportConfiguration(self):
        log = logging.getLogger('bgen_parser::ReportConfiguration')
        log.info(BuildReportLine("BGEN FILE", self.bgen_filename))
        if self.sample_filename is not None:
            log.info(BuildReportLine("SAMPLE FILE", self.sample_filename))
        if self.meta_filename is not None:
            log.info(BuildReportLine("META FILE", self.meta_filename))


    def load_family_details(self, pheno_covar):
        """Load family data updating the pheno_covar with  family ids found.

        :param pheno_covar: Phenotype/covariate object
        :return: None
        """
        libgwas.timer.report_period("Loading Family Details")
        self.sample_ids = list(self.bgen['samples'])

        artificial_ids = False
        if self.sample_ids[0] == "sample_0":
            artificial_ids = True

        mask_components = [0 for x in self.sample_ids]  # 1s indicate an individual is to be masked out
        # Sample file is only necessary if the user chose not to include them
        # in the bgen file. bgen_reader will generate it's own ids if they
        # don't exist, but we don't care to worry about those.
        if self.sample_filename is not None:
            with open(self.sample_filename) as file:
                header = file.readline()
                format = file.readline()        #
                self.file_index = 0

                sample_index = 0

                for line in file:
                    words = line.strip().split()
                    indid = words[0]

                    if artificial_ids or self.sample_ids[sample_index] == indid:
                        if not DataParser.valid_indid(indid):
                            mask_components[sample_index] = 1
                    sample_index += 1

        sample_index = 0
        for indid in self.sample_ids:
            if mask_components[sample_index] == 0:
                pheno_covar.add_subject(indid, phenotype=pheno_covar.missing_encoding)
            sample_index += 1

        self.ind_mask = numpy.array(mask_components, dtype=numpy.int8)
        self.geno_mask = self.ind_mask.reshape(self.ind_mask.shape[0], 1).repeat(3, axis=1)

        self.ind_count = self.ind_mask.shape[0]
        pheno_covar.freeze_subjects()
        libgwas.timer.report_period("Family Details loaded.")
    def open_bgen(self):
        if self.bgen is None:
            print(f"Samples File path: {self.sample_filename}")
            self.bgen = bgen_reader.read_bgen(self.bgen_filename,
                                            metafile_filepath=None,
                                              samples_filepath=self.sample_filename,
                                              verbose=True)
            self.markers_raw = self.bgen['variants'].compute()
            libgwas.timer.report_period("Marker data loaded: %d variants found " % (len(self.markers_raw)))
        self.bgen_idx = self.bgen_start_idx

    def parse_variant(self, index):
        log = logging.getLogger('bgen_parser::open_bgen')
        print("-->", index, self.markers_raw.shape)
        print("--Parse Variants: %d (%d) " % (index, self.markers_raw.shape[0]))
        if index < self.markers_raw.shape[0]:
            locus = Locus()
            v = self.markers_raw.loc[index]

            locus.chr = v.chrom.lstrip("0")

            if locus.chr.strip() == "" and Parser.default_chromosome != -1:
                locus.chr = Parser.default_chromosome
            locus.pos = v.pos
            locus.alleles = v.allele_ids.split(",")
            locus.rsid = v.rsid
            for component in locus.rsid.split(":"):
                if component[0:3] == 'rs':
                    locus.rsid = component
            locus.cur_idx = index
            locus.nalleles = v.nalleles
            #libgwas.timer.report_period("--%s:%d %s" % (str(locus.chr), locus.pos, locus.rsid))
            return locus
        libgwas.timer.report_period("ParseVariant: %d out of loci to consider" % (index))
        raise StopIteration

    # We'll assume that all files that have been associated with this
    # "dataset" are to be considered.
    def load_genotypes(self):
        """Prepares the files for genotype parsing.
        :return: None
        """
        self.open_bgen()
        missing = None
        locus_count = 0
        libgwas.timer.report_period("Loading Genotypes ")

        # identify individual's missingness if the threshold is set
        if DataParser.ind_miss_tol < 1.0:
            for locus in self:
                if DataParser.boundary.BoundaryCompare(BoundaryCheck.get_valid_chrom(iteration.chr), iteration.pos, iteration.rsid) < 0:
                    self.bgen_start_idx += 1
                if missing is None:
                    missing = numpy.zeros(locus.genotype_data.shape[0])

                missing += ((locus.genotype_data == DataParser.missing_representation))
                locus_count += 1

            max_missing = DataParser.ind_miss_tol * locus_count
            dropped_individuals = 0+(max_missing < missing)
            self.ind_mask = self.ind_mask | dropped_individuals

        valid_individuals = numpy.sum(self.ind_mask==0)
        self.min_likely_hets = float(valid_individuals) * DataParser.min_maf

        self.max_missing_geno = DataParser.snp_miss_tol * float(valid_individuals)
        libgwas.timer.report_period("Genotypes Loaded. Starting Index: %d. Locus Count: %d" % (self.bgen_start_idx, locus_count))

    def get_next_line(self):
        """If we reach the end of the file, we simply open the next, until we \
        run out of archives to process"""
        info = 1.0
        locus = self.parse_variant(self.bgen_idx)
        self.bgen_idx += 1
        return locus, info
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
            iteration.alleles = snpdata.alleles
            nalleles = snpdata.nalleles
            iteration.rsid = snpdata.rsid

            if DataParser.boundary.TestBoundary(BoundaryCheck.get_valid_chrom(iteration.chr), iteration.pos, iteration.rsid):
                iteration.genotype_data = numpy.ma.MaskedArray(self.bgen['genotype'][self.bgen_idx - 1].compute(), self.geno_mask).compressed().reshape(-1, 3)
                libgwas.timer.report_period("-  %d %s:%s - Done"% (self.bgen_idx, iteration.chr, str(iteration.pos)))
                # Assuming that if we are missing the first dose, then we are missing them all
                likely_hets = numpy.sum(iteration.genotype_data[:, 1] > Parser.het_threshold)

                # Skip over things that are likely to be fixed loci
                isvalid = likely_hets > self.min_likely_hets
                if isvalid:
                    # technically, there is a ploidy value in there someplace which should end up as 0 for data without information
                    iteration.missing_genotypes = numpy.sum(iteration.genotype_data, axis=1) < 0.1
                    libgwas.timer.report_period("- missingness identified")
                return isvalid

                return True
        return False


    def filter_genotypes(self, missing):
        genotypes = numpy.ma.MaskedArray(self.bgen['genotype'][self.bgen_idx - 1].compute(),
                                         self.geno_mask).compressed().reshape(-1, 3)

        estimate = None
        maf = None
        additive_estimate = genotypes[:, 1] + 2 * genotypes[:, 2]
        if encoding == impute_parser.Encoding.Dominant:
            estimate = genotypes[:, 1] + genotypes[:, 2]
        elif encoding == impute_parser.Encoding.Additive:
            estimate = additive_estimate
        elif encoding == impute_parser.Encoding.Recessive:
            estimate = genotypes[2]
        elif encoding == impute_parser.Encoding.Genotype:
            estimate = numpy.full_like(genotypes.shape, 2)
            estimate[genotypes[:, 1] > genotypes[:, 0] and genotypes[:, 1] > genotypes[:, 2]] = 1
            estimate[genotypes[:, 0] > genotypes[:, 1] and genotypes[:, 0] > genotypes[:, 2]] = 0
        iteration.non_missing_alc = genotypes.shape[0] * 2
        maf = numpy.mean(additive_estimate) / 2
        iteration.allele_count2 = maf * (iteration.non_missing_alc * 2)
        iteration.effa_freq = maf
        iteration.major_allele, iteration.minor_allele = snpdata.alleles
        if maf > 0.5:
            iteration.min_allele_count = iteration.non_missing_alc - iteration.allele_count2
            iteration.maj_allele_count = iteration.allele_count2
            maf = 1.0 - maf
            mallele = iteration.major_allele
            iteration.maj_allele = iteration.minor_allele
            iteration.minor_allele = mallele
        else:
            iteration.min_allele_count = iteration.allele_count2
            iteration.maj_allele_count = iteration.non_missing_alc - iteration.allele_count2
        iteration._maf = maf
        iteration.genotype_data = numpy.array(estimate)



        plocus = AlleleCounts(genotypes)
        plocus.het_count = numpy.sum(genotypes==1)
        a1 = numpy.sum(genotypes==0) + plocus.het
        a2 = numpy.sum(genotypes==2) + plocus.het

        plocus.set_allele_counts(a1, a2, self.alleles)

    def get_effa_freq(self, genotypes):
        return numpy.mean(genotypes)/2

    def __iter__(self):
        """Reset the file and begin iteration"""

        loc = ParsedLocus(self)
        loc._extract_genotypes = gen_dosage_extraction
        return loc
