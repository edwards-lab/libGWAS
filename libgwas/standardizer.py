
import numpy
import pheno_covar

from exceptions import InvariantVar

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

def get_standardizer():
    global _standardizer
    return _standardizer
def set_standardizer(std):
    global _standardizer
    _standardizer = std

class StandardizedVariable(object):
    """Optional plugin object that can be used to standardize covariate and
       phenotype data.

       Many algorithms require that input be standardized in some way in order
       to work properly, however, rescaling the results is algorithm specific.
       In order to facilitate this situation, application authors can
       write up application specific Standardization objects for use with
       the data parsers.

       """
    def __init__(self, pc):
        #: mask representing missingness (1 indicates missing)
        self.missing = []
        #: number of covars
        self.covar_count = len(pc.covariate_data)
        #: number of phenotypes
        self.pheno_count = len(pc.phenotype_data)
        #: Standardized covariate data
        self.covariates = None
        #: standardized phenotype data
        self.phenotypes = None

        for pheno in pc.phenotype_data:
            missing = pheno == pheno_covar.PhenoCovar.missing_encoding
            for idx in range(0, self.covar_count):
                missing = missing | (pc.covariate_data[idx] == pheno_covar.PhenoCovar.missing_encoding)
            self.missing.append(missing)
        #: index of the current phenotype
        self.idx = 0
        #: Reference back to the pheno_covar object for access to raw data
        self.datasource = pc


    def get_variables(self, missing_in_geno=None):
        """Extract the complete set of data based on missingness over all
        for the current locus.

        :param missing_in_geno: mask associated with missingness in genotype
        :return: (phenotypes, covariates, nonmissing used for this set of vars)
        """
        count = 0
        mismatch = 0

        if missing_in_geno is None:
            nonmissing = numpy.invert(self.missing[self.idx])
        else:
            nonmissing = numpy.invert(self.missing[self.idx] | missing_in_geno)
        nmcount = sum(nonmissing)
        covars = numpy.zeros((self.covar_count, nmcount))
        for idx in range(0, self.covar_count):
            covars[idx] = self.covariates[idx][nonmissing]
            min = covars[idx][covars[idx] != pheno_covar.PhenoCovar.missing_encoding].min()
            max = covars[idx][covars[idx] != pheno_covar.PhenoCovar.missing_encoding].max()
            if min == max:
                raise InvariantVar("Covar %s doesn't have enough variation to continue" % (self.datasource.covariate_labels[idx]))

        min = self.phenotypes[self.idx][nonmissing].min()
        max = self.phenotypes[self.idx][nonmissing].max()
        if min == max:
            raise InvariantVar("Phenotype %s doesn't have enough variation to continue" % (self.datasource.phenotype_names[self.idx]))
        return (self.phenotypes[self.idx][nonmissing], covars, nonmissing)

    def get_phenotype_name(self):
        """Returns current phenotype name"""
        return self.datasource.phenotype_names[self.idx]

    def get_covariate_name(self, idx):
        """Return label for a specific covariate

        :param idx: which covariate?
        :return: string label
        """
        return self.datasource.covariate_labels[idx]

    def get_covariate_names(self):
        """Return all covariate labels as a list

        :return: list of covariate names
        """
        return self.datasource.covariate_labels

    def standardize(self):
        """Stub for the appropriate standardizer function

        Each Standardizer object will do it's own thing here.
        """
        pass

    def destandardize(self):
        """Stub for the appropriate destandardizer function.

        Each object type will do it's own thing here.
        """
        pass


class NoStandardization(StandardizedVariable):
    """This is mostly a placeholder for standardizers. Each application will
    probably have a specific approach to standardizing/destandardizing the
    input/output.

    """

    def __init__(self, pc):
        super(NoStandardization, self).__init__(pc)
    def standardize(self):
        """Standardize the variables within a range [-1.0 and 1.0]

        This replaces the local copies of this data. When it's time to
        scale back, use destandardize from the datasource for that.

        """
        self.covariates = self.datasource.covariate_data
        self.phenotypes = self.datasource.phenotype_data

    def destandardize(self, estimates, se, **kwargs):
        """When the pheno/covar data has been standardized, this can be
        used to rescale the betas back to a meaningful value using the
        original data.

        For the "Un-standardized" data, we do no conversion.

        """

        return estimates, se, kwargs["pvalues"]

#: This should be set to an appropriate object by the application
_standardizer = NoStandardization
