import collections
from exceptions import InvalidBoundarySpec
from . import BuildReportLine
import os

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


class BoundaryCheck(object):
    """Record boundary specifications from user to control traversal.

    Default boundaries are specified in numerical positions along a
    single chromosome. Users are permitted to provide boundaries in
    3 forms: Bases, Kilobasees and Megabases. All are recorded as
    single base offsets from the beginning of the chromosome (starting
    at 1).

    The valid setting doesn't mean the boundary object is invalid, only
    that no actual boundary ranges have been provided. This is done to
    allow the user interface code to be a little simpler (i.e. if the
    user didn't provide bounds using numerical boundaries, it can try
    instantiating a SnpBoundary and pass the relevant arguments to
    that object. If none are valid, then either can be used, at which
    point both act as chromosome boundaries or simple SNP filters)

    If chrom is specfied, all SNPs and boundaries are expected to
    reside on that chromosome.

    """
    chrom = -1
    def __init__(self, bp=(None, None), kb=(None, None), mb=(None, None)):
        """Initialize boundary

        :param bp: limit range in base pairs
        :param kb: limit range in kilobases
        :param mb: limit range in megabases
        :return: None

        If any of the range objects contains a valid pair, valid is set to
        be True. Only one boundary pair is accepted. So, if multiples are
        provided, the most specific is accepted (bp is most specific).

        """
        #: List of RS Numbers to be ignored
        self.ignored_rs = []

        #: List of RS Numbers to be targeted (ignors all but those listed)
        self.target_rs = []

        #: Indices of loci that are to be dropped
        #: {chr=>[pos1, pos2, ..., posN]}
        self.dropped_snps   = collections.defaultdict(set)

        #: Actual boundary details in BP
        self.bounds = []

        #: True if boundary conditions remain true
        self.valid = True

        if bp[0] != None and bp[0] + bp[1] > 0:
            self.bounds = bp
        elif kb[0] != None and kb[0] + kb[1] > 0:
            self.bounds = [1000*i for i in kb]
        elif mb[0] != None and mb[0] + mb[1] > 0:
            self.bounds = [1000000*i for i in mb]
        else:
            self.valid = False

        if len(self.bounds) > 0:
            if BoundaryCheck.chrom == -1:
                raise InvalidBoundarySpec(("--chr must be present for " +
                                           "positional filtering to work"))
            # If there is a meaningful boundary configuration but a meaningless
            # chromosome in place, then there is a problem
            try:
                chr = int(BoundaryCheck.chrom)
            except:
                self.valid = False
        #: Is set once the upper limit has been exceeded
        self.beyond_upper_bound = False

    def LoadExclusions(self, snps):
        """ Load locus exclusions.

        :param snps: Can either be a list of rsids or a file containing rsids.
        :return: None

        If snps is a file, the file must only contain RSIDs separated
        by whitespace (tabs, spaces and return characters).
        """

        snp_names = []
        if len(snps) == 1 and os.path.isfile(snps[0]):
            snp_names = open(snps).read().strip().split()
        else:
            snp_names = snps
        for snp in snp_names:
            if len(snp.strip()) > 0:
                self.ignored_rs.append(snp)

    def TestBoundary(self, chr, pos, rsid):
        """Test if locus is within the boundaries and not to be ignored.

        :param chr: Chromosome of locus
        :param pos: BP position of locus
        :param rsid: RSID (used to check for exclusions)
        :return: True if locus isn't to be ignored
        """

        # We can skip over anything that has a negative boundary
        if pos < 0:
            return False

        # Skip over anything that has been explicitly ignored
        if rsid in self.ignored_rs:
            return False

        # If Chromosome isn't defined, then we have no bounds
        if BoundaryCheck.chrom == -1:
            return pos not in self.dropped_snps[chr]

        self.beyond_upper_bound = chr > BoundaryCheck.chrom

        if chr == BoundaryCheck.chrom:
            if len(self.bounds) == 0:
                return pos not in self.dropped_snps[chr]

            if (pos>=self.bounds[0]) and (pos<=self.bounds[1]):
                return pos not in self.dropped_snps[chr]
            self.beyond_upper_bound = pos > self.bounds[1]

        return rsid in self.target_rs

    def NoExclusions(self):
        """Determine that there are no exclusion criterion in play

        :return: True if there is no real boundary specification of any kind.

        Simple method allowing parsers to short circuit the determination of
        missingness, which can be moderately compute intensive.
        """

        if len(self.ignored_rs) + len(self.target_rs) + len(self.bounds) == 0:
            return BoundaryCheck.chrom == -1
        return False

    def ReportConfiguration(self, f):
        """Report the boundary configuration details

        :param f: File (or standard out/err)
        :return: None
        """

        if BoundaryCheck.chrom != -1:
            print >> f, BuildReportLine("CHROM", BoundaryCheck.chrom)
            if len(self.bounds) > 0:
                print >> f, BuildReportLine("SNP BOUNDARY", "-".join(
                    [str(x) for x in self.bounds]))
        if len(self.ignored_rs) > 0:
            print >> f, BuildReportLine("IGNORED RS", ",".join(self.ignored_rs))
        if len(self.target_rs) > 0:
            print >> f, BuildReportLine("TARGET RS", ",".join(self.target_rs))

    def LoadSNPs(self, snps=[]):
        """Define the SNP inclusions (by RSID). This overrides true boundary \
            definition.

        :param snps: array of RSIDs
        :return: None

        This doesn't define RSID ranges, so it throws InvalidBoundarySpec if it
        encounters what appears to be a range (SNP contains a "-")
        """

        for snp in snps:
            bounds = snp.split("-")
            if len(bounds) == 1:
                if bounds[0] != "":
                    self.target_rs.append(bounds[0])
            else:
                raise InvalidBoundarySpec(snp)


