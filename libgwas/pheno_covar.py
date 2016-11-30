import numpy
from exceptions import MalformedInputFile
from exceptions import InvalidSelection
from exceptions import InvariantVar
from exceptions import NoMatchedPhenoCovars
from standardizer import get_standardizer
import sys

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

class PhenoCovar(object):
    """Store both phenotype and covariate data in a single object.

    Provide iterable interface to allow evaluation of multiple phenotypes
    easily. Covariates do not change during iteration. Missing is updated
    according to the missing content within the phenotype (and covariates
    as well).
    """

    #: Do we use sex as a covariate?
    sex_as_covariate = False
    #: Internal encoding for missingness
    missing_encoding = -9

    # Prime with pedigree data and the --sex == True. This will optionally add/activate the first covariate, SEX
    # Load phenotype data from file. If this happens, we'll overwrite the pedigree based data
    # Load covariates from file. This will not replace the sex values pulled from the pedigree file.

    def __init__(self):
        #: Raw phenotype data with every possible phenotype [[ph1],[ph2],etc]
        self.phenotype_data = [[]]
        #: All covariate data [[cov1],[cov2],etc]
        self.covariate_data = []
        #: Pedigree information {FAMID:INDID => index, etc}
        self.pedigree_data = {}
        #: True indicates an individual is to be excluded
        self.individual_mask = []
        #: List of covariate names from header, if provided
        #: SEX is implied, if sex_as_covariate is true.
        #: Covariates loaded without header are simply named Cov-N
        self.covariate_labels = []
        #: List of phenotype names from header, if provided.
        #: If no header is found, the phenotype is simply named Pheno-N
        self.phenotype_names = ["Pheno-1"]
        #: Allows you to turn off standardization
        self.do_standardize_variables = False

        if PhenoCovar.sex_as_covariate:
            self.covariate_labels.append("SEX")
            self.covariate_data.append([])

        #: finalized data ready for analysis
        self.test_variables = None

    def __iter__(self):
        if self.test_variables is None:
            self.prep_testvars()
        if len(self.phenotype_data) == 0 or len(self.phenotype_data[0]) == 0:
            raise StopIteration
        for idx in range(0, len(self.phenotype_data)):
            self.test_variables.idx = idx
            yield self.test_variables
        raise StopIteration

    def prep_testvars(self):
        """Make sure that the data is in the right form and standardized as
        expected.
        """
        self.phenotype_data = numpy.array(self.phenotype_data)
        self.covariate_data = numpy.array(self.covariate_data)
        self.test_variables = get_standardizer()(self)
        self.test_variables.standardize()

    def destandardize_variables(self, tv, blin, bvar, errBeta, nonmissing):
        """Destandardize betas and other components."""
        return self.test_variables.destandardize(tv, blin, bvar, errBeta, nonmissing)


    def freeze_subjects(self):
        """Converts variable data into numpy arrays.

        This is required after all subjects have been added via the
        add_subject function, since we don't know ahead of time who is
        participating in the analysis due to various filtering possibilities.
        """
        self.phenotype_data = numpy.array(self.phenotype_data)
        self.covariate_data = numpy.array(self.covariate_data)



    def add_subject(self, ind_id, sex=None, phenotype=None):
        """Add new subject to study, with optional sex and phenotype

        Throws MalformedInputFile if sex is can't be converted to int
        """

        self.pedigree_data[ind_id] = len(self.phenotype_data[0])
        if phenotype != None:
            if type(self.phenotype_data) is list:
                self.phenotype_data[0].append(phenotype)
            else:
                self.phenotype_data[-1, len(self.individual_mask)] = phenotype
        self.individual_mask.append(0)

        if PhenoCovar.sex_as_covariate:
            try:
                self.covariate_data[0].append(float(sex))
            except Exception, e:
                raise MalformedInputFile("Invalid setting, %s, for sex in pedigree" % (sex))

    def load_phenofile(self, file, indices=[], names=[], sample_file=False):
        """Load phenotype data from phenotype file

        Whitespace delimited, FAMID, INDID, VAR1, [VAR2], etc

        Users can specify phenotypes of interest via indices and names.
        Indices are 1 based and start with the first variable. names
        must match name specified in the header (case is ignored). """

        file.seek(0)
        self.phenotype_names = []
        if file:
            header = ["", "", ""]
            header = file.readline().strip().upper().split()
            line_number = 0
            valid_indices = [int(x) for x in indices]
            valid_names  = []

            for name in names:
                if name.strip() != "":
                    valid_names.append(name)

            # Rows ignored because family/individual missing from our pedigree data
            ignored_data = []
            # We can accept a default phenotype column if we only have 3 columns
            if len(header) == 3:
                if len(valid_names) + len(valid_indices) == 0:
                    valid_indices.append(1)
            if header[0].upper() == "FID":
                phenotype_names = header[2:]
                for name in valid_names:
                    try:
                        valid_indices.append(phenotype_names.index(name.upper()) + 1)
                    except:
                        raise InvalidSelection(
                            "The name, %s, was not found in %s" % (name, file.name)
                        )
                line_number = 1

                if len(valid_indices) > 0 and max(valid_indices) > (len(header)-2):
                    raise InvalidSelection(
                        "The index, %s, is larger than the number of entries in the file, %s:%s" %
                        (max(valid_indices), file.name, line_number)
                    )

                for i in xrange(0, len(valid_indices)):
                    self.phenotype_names.append(phenotype_names[valid_indices[i]-1])
                # Dump the second line for sample_file==True
                if sample_file:
                    file.readline()
                    line_number += 1
            # if we don't have a header, we'll create dummy names
            else:
                file.seek(0)
                if len(valid_names) > 0:
                    raise MalformedInputFile(
                        "Names only work with phenotype files with headers: %s for file %s" %
                        (",".join(names), file.name)
                    )
                if len(valid_indices) > 0 and max(valid_indices) > (len(header)-2):
                    raise InvalidSelection(
                        "The index, %s, is larger than the number of entries in the file, %s:%s" %
                        (max(valid_indices), file.name, line_number)
                    )

                self.phenotype_names = []
                for i in xrange(0, len(valid_indices)):
                    self.phenotype_names.append("Pheno-%s" % (valid_indices[i]))

            pheno_count = len(valid_indices)
            self.phenotype_data = numpy.empty((pheno_count, len(self.pedigree_data)))
            self.phenotype_data.fill(PhenoCovar.missing_encoding)

            for line in file:
                line_number += 1
                words   = line.split()
                iid = ":".join(words[0:2])

                # Indexes are 1 based...silly humans
                if len(valid_indices) > 0 and max(valid_indices) > (len(words)-2):
                    raise InvalidSelection(
                        "The index, %s, is larger than the number of entries in the file, %s:%s" %
                        (max(valid_indices), file.name, line_number)
                    )


                pidx = 0
                for idx in valid_indices:
                    try:
                        pheno = float(words[1+idx])
                        if iid not in self.pedigree_data:
                            ignored_data.append(iid)
                        else:
                            self.phenotype_data[pidx][self.pedigree_data[iid]] = pheno
                        pidx += 1
                    except:
                        raise MalformedInputFile(
                            ("Invalid input found in phenotype file on line: %s:%d. \n"+
                            "The line in question looks like this: \n--> %s") %
                            (file.name, line_number, line.strip())
                        )
        if self.phenotype_data.shape[1] == len(ignored_data):
            raise NoMatchedPhenoCovars("No matching individuals were found in the phenotype file")

    def load_covarfile(self, file, indices=[], names=[], sample_file=False):
        """Load covariate data from file.

        Unlike phenofiles, if we already have data, we keep it (that would be the sex covariate)"""

        # Clean up input in case we are given some empty values
        var_indices = []
        for x in indices:
            try:
                var_indices.append(int(x) + 1)
            except:
                pass

        var_names = []
        for name in names:
            if name.strip() != "":
                var_names.append(name)

        if PhenoCovar.sex_as_covariate:
            self.covariate_labels = ["SEX"]
        file.seek(0)
        if file:
            header = file.readline().strip().split()
            line_number = 0

            # Rows ignored because family/individual missing from our pedigree data
            ignored_data = []
            # We can accept a default phenotype column if we only have 3 columns
            if len(header) == 3:
                if len(var_names) + len(var_indices) == 0:
                    var_indices.append(2)

            if header[0].upper() == "FID":
                for name in var_names:
                    if name.strip() != "":
                        try:
                            var_indices.append(header.index(name))
                        except:
                            raise InvalidSelection("The name, %s, was not found in %s" % (name, file.name))
                line_number = 1


                if len(var_indices) > 0 and max(var_indices) > (len(header)):
                    raise InvalidSelection("The index, %s, is larger than the number of entries in the file, %s:%s" % (max(var_indices), file.name, line_number))

                for index in var_indices:
                    self.covariate_labels.append(header[index])
                # Dump the second line for sample_file==True
                if sample_file:
                    file.readline()
            else:
                file.seek(0)
                if len(var_names) > 0:
                    raise MalformedInputFile("Names only work with covariate files with headers: %s" % (",".join(names), file.name))
                if len(var_indices) > 0 and max(var_indices) > (len(header)):
                    raise InvalidSelection("The index, %s, is larger than the number of entries in the file, %s:%s" % (max(var_indices), file.name, line_number))

                for i in xrange(0, len(var_indices)):
                    self.covariate_labels.append("Cov-%s" % (var_indices[i]-1))

            covar_data = numpy.empty((len(var_indices)+len(self.covariate_data), len(self.pedigree_data)))
            covar_data.fill(PhenoCovar.missing_encoding)

            # We have to be careful to keep the sex covariate data if it came from the pedigree
            if PhenoCovar.sex_as_covariate:
                covar_data[0] = self.covariate_data[0]
            self.covariate_data = covar_data

            for line in file:
                line_number += 1
                words   = line.split()
                iid = ":".join(words[0:2])

                # Indexes are 1 based...silly humans
                if len(var_indices) > 0 and max(var_indices) > (len(words)):
                    raise InvalidSelection("The index, %s, is larger than the number of entries in the file, %s:%s" % (max(var_indices), file.name, line_number))

                cidx = 0
                if PhenoCovar.sex_as_covariate:
                    cidx += 1
                for idx in var_indices:
                    cov = float(words[idx])
                    if iid not in self.pedigree_data:
                        ignored_data.append(iid)
                    else:
                        self.covariate_data[cidx, self.pedigree_data[iid]] = cov
                    cidx += 1


        for covar in self.covariate_data:
            if sum(covar !=PhenoCovar.missing_encoding ) == 0:
                raise NoMatchedPhenoCovars("No matching individuals were found in the covariate file")
