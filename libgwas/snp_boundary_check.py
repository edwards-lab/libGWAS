import collections
from exceptions import InvalidBoundarySpec
from . import BuildReportLine
import os
from boundary import BoundaryCheck

__copyright__ = "Eric Torstenson"
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

class SnpBoundaryCheck(BoundaryCheck):
    """RS (or other name) based boundary checking.

    Same rules apply as those for BoundaryCheck, except users can
    provide multiple RS boundary regions. Though, all boundary groups
    must reside on a single chromosome.

    Class members (these are not intended for public consumption):
        * start_bounds    bp location for boundary starts \
                            Currently, only one boundary is permitted. This \
                            is to remain consistant with plink
        * end_bounds      bp location for boundary end (inclusive)
        * ignored_rs      List of RS numbers to be ignored
        * target_rs       List of RS numbers to be targeted
        * dropped_snps    indices of loci that are to be dropped
                            {chr=>[pos1, pos2, ...]}
        * end_rs          This is used during iteration to identify when to
                            turn "off" the current boundary group

    """

    def __init__(self, snps=[]):
        """Initializes everything as empty, and processes snps

        :param snps: List of rsids to add to the target list

        snps may contain singular locus names or boundary pairs
        separated by a dash "-". More than one SNP must be separated by
        comma.

        """
        self.start_bounds   = []
        self.end_bounds     = []
        self.target_rs      = []
        self.ignored_rs     = []
        self.beyond_upper_bound = False

        # We'll drop SNPs by chromosome: chr=>[pos, pos]
        self.dropped_snps   = collections.defaultdict(set)
        self.end_rs         = None  # If this is non-null, we are within a bounded set
        for snp in snps:
            if snp.strip() != "":
                bounds = snp.split("-")
                bound_size = len(bounds)
                if bound_size == 1:
                    self.target_rs.append(bounds[0])
                elif bound_size==2:
                    self.start_bounds.append(bounds[0])
                    self.end_bounds.append(bounds[1])
                else:
                    print "Invalid Bound Size of %s at locus: %s" % (bound_size, snp)
                    raise InvalidBoundarySpec(snp)
        self.valid          = len(self.target_rs)+len(self.start_bounds) > 0
    def ReportConfiguration(self, f):
        """Report the boundary configuration details

        :param f: File (or standard out/err)
        :return: None
        """

        if BoundaryCheck.chrom != -1:
            print >> f, BuildReportLine("CHROM", BoundaryCheck.chrom)
            if len(self.start_bounds) > 0:
                bounds = ",".join(["%s-%s" % (a[0], a[1]) for a in zip(self.start_bounds, self.end_bounds)])
                print >> f, BuildReportLine("SNP BOUNDARY", bounds)
        if len(self.ignored_rs) > 0:
            print >> f, BuildReportLine("IGNORED RS", ",".join(self.ignored_rs))
        if len(self.target_rs) > 0:
            print >> f, BuildReportLine("TARGET RS", ",".join(self.target_rs))



    def TestBoundary(self, chr, pos, rsid):
        """Test if locus is within the boundaries and not to be ignored.

        :param chr: Chromosome of locus
        :param pos: BP position of locus
        :param rsid: RSID (used to check for exclusions)
        :return: True if locus isn't to be ignored
        """

        if self.chrom != -1 and self.chrom != chr:
            return False

        if self.end_rs:
            valid = not (rsid in self.ignored_rs or pos in self.dropped_snps[chr])

            # Be sure to clear our end if we have reached it
            if rsid==self.end_rs:
                self.end_rs = None
            return valid
        else:
            if rsid in self.start_bounds:
                try:
                    idx = self.start_bounds.index(rsid)
                    self.end_rs = self.end_bounds[idx]
                    return rsid not in self.ignored_rs and pos not in self.dropped_snps[chr]

                except:
                    pass
        if rsid in self.target_rs:
            return True
        return False

    def NoExclusions(self):
        """Determine that there are no exclusion criterion in play

        :return: True if there is no real boundary specification of any kind.

        Simple method allowing parsers to short circuit the determination of
        missingness, which can be moderately compute intensive.
        """
        if len(self.start_bounds) + len(self.target_rs) + len(self.ignored_rs) == 0:
            return BoundaryCheck.chrom == -1
        return False
