from locus import Locus

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

class ParsedLocus(Locus):
    """Locus data representing current iteration from a dataset

    Provide an iterator interface for all dataset types.

    """

    def __init__(self, datasource, index=-1):
        """Basic initialization (nothing is currently valid)"""

        super(ParsedLocus, self).__init__()
        #: Reference back to the parser that generated this object
        self.__datasource       = datasource
        #: Index within the list of loci being analyzed
        self.cur_idx            = index
        #: Actual genotype data for this locus
        self.genotype_data      = None

    def next(self):
        """Move to the next valid locus.

        Will only return valid loci or exit via StopIteration exception

        """
        while True:
            self.cur_idx += 1
            if self.__datasource.populate_iteration(self):
                return self

        raise StopIteration


    def __iter__(self):
        """Basic iterator functionality. """

        return self
