

from parsed_locus import ParsedLocus
from locus import Locus

from boundary import BoundaryCheck
import numpy


__copyright__ = "Todd Edwards, Chun Li & Eric Torstenson"
__license__ = "GPL3.0"
#     This file is part of pyGWAS.
#
#     pyGWAS is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     pyGWAS is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with MVtest.  If not, see <http://www.gnu.org/licenses/>.

def check_inclusions(item, included=[], excluded=[]):
    """Everything passes if both are empty, otherwise, we have to check if \
    empty or is present."""
    if (len(included) == 0):
        if len(excluded) == 0 or item not in excluded:
            return True
        else:
            return False
    else:
        if item in included:
            return True
    return False


class DataParser(object):
    """Abstract representation of all dataset parsers

    """

    #: this can be used to filter out loci with too few minor alleles
    min_maf         = 0.00

    #: filter out if a minor allele frequency exceeds this value
    max_maf         = 1.00

    #: Filter SNPs with too many missing
    snp_miss_tol    = 1.0

    #: Filter individuals with too many missing
    ind_miss_tol    = 1.0

    #: Boundary object specifying valid region for analysis
    boundary        = BoundaryCheck()

    #: Filter out specific individuals by individual ID
    ind_exclusions  = []

    #: Filter in specific individuals by individual ID
    ind_inclusions  = []

    #: When false, pedigree header expects no sex column
    has_sex         = True

    #: When false, pedigree header expects no parents columns
    has_parents     = True

    #: When false, pedigree header expects no family id column
    has_fid         = True

    #: When false, pedigree header expects no phenotype column
    has_pheno       = True

    #: When false, pedigree header expects no liability column
    has_liability   = False

    #: External representation of missingness
    missing_representation      = '0'

    # Internal representation of missingness
    missing_storage             = -1

    #: When true, assume that standard pedigree and transposed pedigree are
    #: compressed with gzip
    compressed_pedigree = False

    def get_effa_freq(self, genotypes):
        return numpy.sum(numpy.array(genotypes)-1)/len(genotypes)

    def __iter__(self):
        """Iteration is performed by ParsedLocus"""
        if DataParser.boundary.beyond_upper_bound:
            raise StopIteration

        return ParsedLocus(self)

    @staticmethod
    def valid_indid(indid):
        return check_inclusions(indid, DataParser.ind_inclusions,
                                DataParser.ind_exclusions)


    def get_loci(self):
        loci = []
        for locus in self:
            loci.append(Locus(locus))

            locus.genotype_data = None

        return loci
