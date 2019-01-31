import numpy
import data_parser

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

class AlleleCounts(object):
    """Counter associated with alleles with consideration for missingness
       across all components (genotypes, phenotypes and covariates)
       This will be used for filtering and will provide a rudimentary
       wrapper around a test to ensure that phenotypes are only analyzed
       in the event that the counts pass muster under the user's
       restrictions.

        -- genotypes will be the filtered genotype array which will be used
           for analysis
    """

    def __init__(self, genotypes, alleles, non_missing):
        self.genotypes = genotypes
        self.effa_freq = None
        self.a1_count = 0
        self.a2_count = 0
        self.alleles = alleles
        self.minor_allele = alleles[1]
        self.major_allele = alleles[0]
        self.total_alleles = 0
        self.maf = 0
        self.effa_freq = 0

        self.het_count = 0
        self.missing = numpy.sum(non_missing==0)

    @property
    def hetero_freq(self):
        return float(self.het_count) / (self.total_alleles/2)

    @property
    def freq_missing(self):
        return self.missing / (self.total_alleles/2)

    def set_allele_counts(self, a1, a2, het):
        self.a1_count = a1
        self.a2_count = a2
        self.het_count = het
        self.total_alleles = (a1 + a2)

        self.effa_freq = a2 / float(self.total_alleles)

        if self.effa_freq > 0.5:
            self.maf = a1 / float(self.total_alleles)
            self.minor_allele, self.major_allele = self.alleles
        else:
            self.maf = self.effa_freq

    def __str__(self):
        return " ".join([str(x) for x in [self.a1_count, self.a2_count, self.total_alleles, self.effa_freq, self.maf, self.minor_allele]])
