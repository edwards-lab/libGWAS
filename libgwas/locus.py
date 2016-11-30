
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

import sys
import traceback

class Locus(object):
    def __init__(self, other=None):
        if other:
            #: Chromosome
            self.chr = other.chr
            #: BP Position
            self.pos = other.pos
            #: RSID
            self.rsid = other.rsid
            #: List of alleles present
            self.alleles = list(other.alleles)
            #: total count of heterozygotes observed
            self.hetero_count = other.hetero_count
            #: total number of minor alleles observed
            self.min_allele_count = other.min_allele_count
            #: total number of major alleles observed
            self.maj_allele_count = other.maj_allele_count
            #: total number of missing alleles were observed
            self.missing_allele_count = other.missing_allele_count
            #: minor allele if allele info isn't available
            self._maf = other._maf
            self.cur_idx = other.cur_idx
        else:
            #: minor alelel if allele info isn't available
            self._maf = None
            self.chr = -1
            self.pos = -1
            self.rsid = ""
            self.alleles = ["",""]
            self.hetero_count = -1
            self.min_allele_count = -1
            self.maj_allele_count = -1
            self.missing_allele_count = -1
            self.cur_idx = -1

    def flip(self):
        """This will switch major/minor around, regardless of frequency truth.

        This is intended for forcing one of two populations to relate correctly
        to the same genotype definitions. When flipped, Ps and Qs will be
        backward, and the maf will no longer relate to the "minor" allele
        frequency. However, it does allow clients to use the same calls for each
        population without having to perform checks during those calculations.
        """

        maj_count = self.maj_allele_count
        self.maj_allele_count = self.min_allele_count
        self.min_allele_count = maj_count

        alleles = self.alleles
        self.alleles = [alleles[1], alleles[0]]

    @property
    def sample_size(self):
        """Returns to total sample size"""
        return self.total_allele_count * 0.5

    @property
    def total_allele_count(self):
        """Returns the total number of alleles"""
        return self.min_allele_count + self.maj_allele_count

    @property
    def hetero_freq(self):
        """Returns the frequency of observed heterozygotes (not available with \
        all parsers)"""
        return self.hetero_count / (0.5 * (self.min_allele_count +
                                           self.maj_allele_count))
    @property
    def exp_hetero_freq(self):
        """Returns the estimated frequency of heterozygotes"""
        return 2 * self.p * self.q

    @property
    def p(self):
        """Frequency for first allele"""
        return self.maj_allele_count / float(self.total_allele_count)

    @property
    def q(self):
        """Frequency for second allele"""
        return self.min_allele_count / float(self.total_allele_count)

    @property
    def major_allele(self):
        """Sets/Returns the encoding for the major allele (A, C, G, T, etc)"""
        return self.alleles[0]
    @major_allele.setter
    def major_allele(self, allele):
        self.alleles[0] = allele


    @property
    def minor_allele(self):
        """Sets/Returns the encoding for minor allele"""
        return self.alleles[1]
    @minor_allele.setter
    def minor_allele(self, allele):
        self.alleles[1] = allele

    @property
    def maf(self):
        """Returns the MAF. This is valid for all parsers"""
        if self._maf:
            return self._maf
        return self.q

    def __hash__(self):
        return hash((self.chr, self.pos))

    def __eq__(self, other):
        return (self.chr, self.pos) == (other.chr, other.pos)

    def __neq__(self, other):
        return not(self == other)

    def __cmp__(self, other):
        if self.chr == other.chr:
            return self.pos.__cmp__(other.pos)
        return self.chr.__cmp__(other.chr)

    def __str__(self):
        return "%d\t%d:%d %s %s %s %0.4f %0.4f" % (self.cur_idx, self.chr,
                                                          self.pos, self.rsid,
                                                          self.major_allele,
                                                          self.minor_allele,
                                                          self.maf,
                                                          self.hetero_freq,
-                                                         self.genotype_data)
