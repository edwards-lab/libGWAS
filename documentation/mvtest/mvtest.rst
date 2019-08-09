What is MVtest?
===============

*TODO: Write some background information about the application and it's
scientific basis.*

Documentation
+++++++++++++
Documentation for MVtest is still under construction. However, the application
provides reasonable inline help using standard unix help arguments:

> `mvtest.py -h`

or

> `mvtest.py --help`

In general, overlapping functionality should mimic that of PLINK.

Command-Line Arguments
++++++++++++++++++++++
Command line arguments used by MVtest often mimick those used by PLINK, except
where there is no matching functionality (or the functionality differs
significantly.)

For the parameters listed below, when a parameter requires a value, the value
must follow the argument with a single space separating the two (no '=' signs.)
For flags with no specified value, passing the flag indicates that condition
is to be "activated".

When there is no value listed in the "Type" column, the arguments are *off* by
default and *on* when the argument is present (i.e. by default, compression
is turned off except when the flag, --compression, has been provided.)

Getting help
------------

.. option:: -h, --help

    Show this help message and exit.

.. option:: -v

    Print version number

.. raw:: latex

    \newpage

Input Data
----------
MVtest attempts to mimic the interface for PLINK where appropriate.

All input files should be whitespace delimited. For text based allelic
annotations, 1|2 and A|C|G|T annotation is sufficient. All data must be
expressed as alleles, not as genotypes (except for IMPUTE output, which is a
specialized format that is very different from the other forms).

For Pedigree, Transposed Pedigree and PLINK binary pedigree files, the using
the PREFIX arguments is sufficient and recommended if your files follow the
standard naming conventions.

Pedigree Data
^^^^^^^^^^^^^
Pedigree data is fully supported, however it is not recommended. When loading
pedigree data, MVtest must load the entire dataset into memory prior to
analysis, which can result in a substantial amount of memory overhead that is
unnecessary.

Flags like --no-pheno and --no-sex can be used in any combination creating
MAP files with highly flexible header structures.

.. option:: --file <prefix>

    (filename prefix)
    Prefix for .ped and .map files

.. option:: --ped <filename>

    PLINK compatible .ped file

.. option:: --map <filename>

    PLink compatible .map file

.. option:: --map3

    Map file has only 3 columns

.. option:: --no-sex

    Pedigree file doesn't have column 5 (sex)

.. option:: --no-parents

    Pedigree file doesn't have columns 3 and 4 (parents)

.. option:: --no-fid

    Pedgiree file doesn't have column 1 (family ID)

.. option:: --no-pheno

    Pedigree file doesn't have column 6 (phenotype)

.. option:: --liability

    Pedigree file has column 7 (liability)

.. raw:: latex

    \newpage

PLINK Binary Pedigree
^^^^^^^^^^^^^^^^^^^^^
This format represents the most efficient storage for large GWAS datasets,
and can be used directly by MVtest. In addition to a minimal overhead, plink
style bed files will also run very quickly, due to the efficient disk layout.

.. option:: --bfile <prefix>

    (filename prefix)
    <prefix> for .bed, .bim and .fam files

.. option:: --bed <filename>

    Binary Ped file(.bed)

.. option:: --bim <filename>

    Binary Ped marker file (.bim)

.. option:: --fam <filename>

    Binary Ped family file (.fam)



Transposed Pedigree Data
^^^^^^^^^^^^^^^^^^^^^^^^
Transposed Pedigree data is similar to standard pedigree except that the data
is arranged such that the data is organized as SNPs as rows, instead of
individuals. This allows MVtest to run it's analysis without loading the
entire dataset into memory.

.. option:: --tfile <prefix>

    Prefix for .tped and .tfam files

.. option:: --tped <filename>

    Transposed Pedigre file (.tped)

.. option:: --tfam <filename>

    Transposed Pedigree Family file (.tfam)

Pedigree/Transposed Pedigree Common Flags
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
By default, Pedigree and Transposed Pedigree data is assumed to be uncompressed.
However, MVtest can directly use gzipped data files if they have the extension
.tgz with the addition of the --compressed argument.

.. option:: --compressed

    Indicate that ped/tped files have been compressed with gzip and are named
    with extensions such as .ped.tgz and .tped.tgz

.. raw:: latex

    \newpage

IMPUTE output
^^^^^^^^^^^^^
MVtest doesn't call genotypes when performing analysis and allows users to
define which model to use when analyzing the data. Due to the fact that there
is no specific location for chromosome within the input files, MVtest requires
that users provide chromosome, impute input file and the corresponding .info
file for each imputed output.

Due to the huge number of expected loci, MVtest allows users to specify an
offset and file count for analysis. This is to allow users to run multiple
jobs simultaneously on a cluster and work individually on separate impute
region files. Users can segment those regions even further using standard
MVtest region selection as well.

By default, all imputed data is assumed to be compressed using gzip.

Default naming convention is for impute data files to end in .gen.gz and
the info files to have the same name except for the end being replaced
by .info.

.. option:: --impute <filename>

    File containing list of impute output for analysis

.. option:: --impute-fam <filename>

    File containing family details for impute data

.. option:: --impute-offset <integer>

    Impute file index (1 based) to begin analysis

.. option:: --impute-count <integer>

    Number of impute files to process (for this node). Defaults to all remaining.

.. option:: --impute-uncompressed

    Indicate that the impute input is not gzipped, but plain text

.. option:: --impute-encoding

    (additive,dominant or recessive)

    Genetic model to be used when analyzing imputed data.

.. option:: --impute-info-ext <extension>

    Portion of filename denotes info filename

.. option:: --impute-gen-ext <extension>

    Portion of filename that denotes gen file

.. option:: --impute-info-thresh <float>

    Threshold for filtering imputed SNPs with poor 'info' values

IMPUTE File Input
^^^^^^^^^^^^^^^^^
When performing an analysis on IMPUTE output, users must provide a single file
which lists each of the gen files to be analyzed. This plain text file contains
2 (or optionally 3) columns for each gen file:


.. tabularcolumns:: |p{5cm}|p{5cm}|p{5cm}|
================  ==============  ===============================
 **Chromosome**    **Gen File**    **.info <filename> (optional)**
================  ==============  ===============================
  N                <filename>         <filename>
  ...              ...              ...
================  ==============  ===============================

The 3rd column is only required if your .info files and .gen files are not
the same except for the <extension>.

.. raw:: latex

    \newpage

MACH output
^^^^^^^^^^^
Users can analyze data imputed with MACH. Because most situations require
many files, the format is a single file which contains either pairs of
dosage/info files, or, if the two files share the same filename except for
extensions, one dosage file per line.

.. important::

    MACH doesn't provide anywhere to store chromosome and positions. Users may
    wish to embed this information into the first column inside the .info file.
    Doing so will allow MVtest to recognize those values and populate the
    corresponding fields in the report.

    To use this feature, users much use the --mach-chrpos field and their ID
    columns inside the .info file must be formatted in the following way:

    chr:pos (optionally :rsid)

    When the --mach-chrpos flag is used, MVtest will fail when it encounters
    IDs that aren't in this format and there must be at least 2 'fields' (i.e.
    there must be at least one ":" character.

    When processing MACH imputed data without this special encoding of IDs,
    MCtest will be unable to recognize positions. As a result, unless the
    --mach-chrpos flag is present, MVtest will exit with an error if the
    user attempts to use positional filters such as --from-bp, --chr, etc. 

When running MVtest using MACH dosage on a cluster, users can instruct a given
job to anlyze data from a portion of the files contained within the MACH
dosage file list by changing the --mach-offset and --mach-count arguments. By
default, the offset starts with 1 (the first file in the dosage list) and runs
all it finds. However, if one were to want to split the jobs up to analyze
three dosage files per job, they might set those values to --mach-offset 1
--mach-count 3 or --mach-offset 4 --mach-count 3 depending on which job
is being defined.

In order to minimize memory requirements, MACH dosage files can be loaded
incrementally such that only N loci are stored in memory at a time. This can
be controlled using the --mach-chunk-size argument. The larger this number is,
the faster MVtest will run (fewer times reading from file) but the more
memory is required.

.. option:: --mach <filename>

    File containing list of dosages, one per line. Optionally, lines may
    contain the info names as well (separated by whitespace) if the two
    <filename>s do not share a common base name.

.. option:: --mach-offset <integer>

    Index into the MACH file to begin analyzing

.. option:: --mach-count <integer>

    Number of dosage files to analyze

.. option:: --mach-uncompressed

    By default, MACH input is expected to be gzip compressed. If data is plain
    text, add this flag. *It should be noted that dosage and info files should
    be either both compressed or both uncompressed.*

.. option:: --mach-chunk-size <integer>

    Due to the individual orientation of the data, large dosage files are parsed
    in chunks in order to minimize excessive memory during loading

.. option:: --mach-info-ext <extension>

    Indicate the <extension> used by the mach info files

.. option:: --mach-dose-ext <extension>

    Indicate the <extension> used by the mach dosage files

.. option:: --mach-min-rsquared <float>

    Indicate the minimum threshold for the rsqured value from the .info files
    required for analysis.

.. option:: --mach-chrpos

    When set, MVtest expects IDs from the .info file to be in the format
    chr:pos:rsid (rsid is optional). This will allow the report to contain
    positional details, otherwise, only the RSID column will have a value
    which will be the contents of the first column from the .info file

MACH File Input
^^^^^^^^^^^^^^^
When running an analysis on MACH output, users must provide a single file which
lists of each dosage file and (optionally) the matching .info file. This file
is a simple text file with either 1 column (the dosage filename) or 2 (dosage
filename followed by the info filename separated by whitespace).

The 2nd column is only required if the filenames aren't identical except for
the extension.

.. tabularcolumns:: |p{6cm}|p{9cm}|
==============================  ==================================
**Col 1 (dosage <filename>)**     **Col 2 (optional info <filename>)**
==============================  ==================================
  <filename>.dose                   <filename>.info
  ...                             ...
==============================  ==================================

Phenotype/Covariate Data
^^^^^^^^^^^^^^^^^^^^^^^^
Phenotypes and Covariate data can be found inside either the standard pedigree
headers or within special PLINK style covariate files. Users can specify
phenotypes and covariates using either header names (if a header exists in
the file) or by 1 based column indices. An index of 1 actually means the
first variable column, not the first column. In general, this will be the
3rd column, since columns 1 and 2 reference FID and IID.

.. option:: --pheno <filename>

    File containing phenotypes. Unless --all-pheno is present, user must
    provide either index(s) or label(s) of the phenotypes to be analyzed.

.. option:: --mphenos LIST

    Column number(s) for phenotype to be analyzed if number of columns > 1.
    Comma separated list if more than one is to be used.

.. option:: --pheno-names LIST

    Name for phenotype(s) to be analyzed (must be in --pheno file). Comma
    separated list if more than one is to be used.

.. option:: --covar <filename>

    File containing covariates

.. option:: --covar-numbers LIST

    Comma-separated list of covariate indices

.. option:: --covar-names LIST

    Comma-separated list of covariate names

.. option:: --sex

    Use sex from the pedigree file as a covariate

.. option:: --missing-phenotype CHAR

    Encoding for missing phenotypes as can be found in the data.

.. option:: --all-pheno

    When present, mv-test will run each phenotypes found inside the phenotype
    file.

.. raw:: latex

    \newpage

Restricting regions for analysis
--------------------------------
When specifying a range of positions for analysis, a chromosome must be present.
If a chromosome is specified but is not accompanied by a range, the entire
chromosome will be used. Only one range can be specified per run.

In general, when specifying region limits, --chr must be defined unless using
generic MACH input (which doesn't define a chromosome number nor position, in
which case positional restrictions do not apply).

.. option:: --snps LIST

    Comma-delimited list of SNP(s): rs1,rs2,rs3-rs6

.. option:: --chr <integer>

    Select Chromosome. If not selected, all chromosomes are to be analyzed.

.. option:: --from-bp <integer>

    SNP range start

.. option:: --to-bp <integer>

    SNP range end

.. option:: --from-kb <integer>

    SNP range start

.. option:: --to-kb <integer>

    SNP range end

.. option:: --from-mb <integer>

    SNP range start

.. option:: --to-mb <integer>

    SNP range end

.. option:: --exclude LIST

    Comma-delimited list of rsids to be excluded

.. option:: --remove LIST

    Comma-delimited list of individuals to be removed from analysis. This must

    be in the form of family_id:individual_id

.. option:: --maf <float>

    Minimum MAF allowed for analysis

.. option:: --max-maf <float>

    MAX MAF allowed for analysis

.. option:: --geno <integer>

    MAX per-SNP missing for analysis

.. option:: --mind <integer>

    MAX per-person missing

.. option:: --verbose

    Output additional data details in final report