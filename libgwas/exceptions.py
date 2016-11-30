
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

class ReportableException(Exception):
    """Simple exeception with message"""
    def __init__(self, msg):
        self.msg = msg


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


class InvalidSelection(MalformedInputFile):
    """Indicate that the user provided input that is meaningless.

    This is likely a situation where the user provided an invalid name
    for a phenotype or covariate. Probably a misspelling.

    """
    pass
