
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


class ReportableException(Exception):
    """Simple exeception with message"""
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg

    def __repr__(self):
        return self.msg

class UnsolvedLocus(ReportableException):
    def __init__(self, msg):
        super(UnsolvedLocus, self).__init__(msg)

class TooManyAlleles(ReportableException):
    """Indicate locus found with more than 2 alleles"""
    def __init__(self, chr=None, rsid=None, pos=None, alleles=None, index=None,
                 prefix="Too many alleles: "):
        #: Chromosome
        self.chr = chr

        #: BP Position
        self.pos = pos

        #: RSID
        self.rsid = rsid

        #: Allele 1 and 2
        self.alleles = alleles

        #: Index of the locus within the file
        self.index = index

        super(TooManyAlleles, self).__init__(
            "%s %s:%s (%s) %s" %
            (prefix, self.chr, self.pos, self.rsid, self.alleles))

class TooMuchMissingpPhenoCovar(ReportableException):
    def __init__(self, pheno_name, missing_pct):

        #: Index of the pheno
        self.name = pheno_name

        #: Missing Percentag
        self.pct = missing_pct

        msg = "The phenotype %s has too much missing (%.4f)" % (self.name, self.pct)
        super(TooMuchMissingpPhenoCovar, self).__init__(msg)


class TooMuchMissing(ReportableException):
    def __init__(self, chr=None, pos=None, rsid=None, maf=-1.0, miss=-1.0):

        #: Chromosome
        self.chr = chr

        #: BP Position
        self.pos = pos

        #: RSID
        self.rsid = rsid

        #: Frequency that triggered issue
        self.maf = maf

        self.miss = miss

        super(TooMuchMissing, self).__init__(
            "The locus %s:%s (%s) miss freq %0.6f is too large" %
            (self.chr, self.pos, self.rsid, self.miss)
        )

class InvalidFrequency(ReportableException):
    def __init__(self, chr=None, pos=None, rsid=None, maf=-1.0, almin=None, almaj=None, a2c=-1):

        #: Chromosome
        self.chr = chr

        #: BP Position
        self.pos = pos

        #: RSID
        self.rsid = rsid

        #: Frequency that triggered issue
        self.maf = maf

        #: minor allele
        self.minor_allele = almin

        #: major allele
        self.major_allele = almaj

        #: allele 2 count
        self.a2_count = a2c

        super(InvalidFrequency, self).__init__(
            "The locus %s:%s (%s) freq %0.6f is invalid" %
            (self.chr, self.pos, self.rsid, self.maf)
        )

class NanInResult(ReportableException):
    """NaN found in result"""
    def __init__(self, msg = ""):
        super(NanInResult, self).__init__(msg)

class NoMatchedPhenoCovars(ReportableException):
    """No ids matched between pheno or covar and the family data"""
    def __init__(self, msg = ""):
        super(NoMatchedPhenoCovars, self).__init__(msg)

class InvariantVar(ReportableException):
    """No minor allele found"""
    def __init__(self, msg=""):
        super(InvariantVar, self).__init__(msg)

class TooFewAlleles(TooManyAlleles):
    """Indicate fixed allele was found"""
    def __init__(self, chr=None, rsid=None, pos=None, alleles=None, index=None):
        super(TooFewAlleles, self).__init__(chr, rsid, pos, alleles, index,
                                            "Too few alleles: ")


class InvalidBoundarySpec(ReportableException):
    """Indicate boundary specification was malformed or non-sensical"""
    def __init__(self, malformed_boundary):
        self.malformed_boundary = malformed_boundary
        super(InvalidBoundarySpec, self).__init__("%s" % (malformed_boundary))


class MalformedInputFile(ReportableException):
    """Error encountered in data from an input file"""
    def __init__(self, msg):
        super(MalformedInputFile, self).__init__(msg)


class InvalidChromosome(ReportableException):
    def __init__(self, chr):
        self.chr = chr
        msg = "An invalid chromosome specification has been encountered: %s" % (str(chr))
        super(InvalidChromosome, self).__init__(msg)

class InvalidSelection(MalformedInputFile):
    """Indicate that the user provided input that is meaningless.

    This is likely a situation where the user provided an invalid name
    for a phenotype or covariate. Probably a misspelling.

    """
    pass
