mvmany Helper script
=====================

In addition to the analysis program, mvtest.py, a helper script, mvmany.py is
also included and can be used to split large jobs into smaller ones suitable
for running on a compute cluster. Users simply run mvmany.py just like they
would run mvtest.py but with a few additional parameters, and mvmany.py will
build multiple job scripts to run the jobs on multiple nodes. It records most
arguments passed to it and will write them to the scripts that are produced.

It is important to note that mvmany.py simply generates cluster scripts and
does not submit them.

The Default Template
++++++++++++++++++++

When mvmany.py is first run, it will generate a copy of the default template
inside the user's home directory named .mv-many.template. This template is
used to define the job details that will be written to each of the job scripts.
By default, the template is configured for the SLURM cluster software, but can
easily be changed to work with any cluster software that works similarly to
the SLURM job manager, such as TORQUE/PBS or sungrid.

In addition to being able to replace the preprocessor definitions to work with
different cluster manager software, the user can also add user specific
definitions, such as email notifications or account specification, giving the
user the the options necessary to run the software under many different system
configurations.

Example Template (SLURM)
^^^^^^^^^^^^^^^^^^^^^^^^

An example template might look like the following ::

    #!/bin/bash
    #SBATCH --job-name=$jobname
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=$memory
    #SBATCH --time=$walltime
    #SBATCH --error $logpath/$jobname.e
    #SBATCH --output $respath/$jobname.txt

    cd $pwd

    $body

It is important to note that this block of text contains a mix of SLURM
preprocessor settings (such as #SBATCH --job-name) as well as variables
which will be replaced with appropriate values (such as $jobname being replaced
with a string of text which is unique to that particular job). Each cluster
type has it's own syntax for setting the necessary variables and it is assumed
that the user will know how to correctly edit the default template to suit
their needs.

Example TORQUE Template
^^^^^^^^^^^^^^^^^^^^^^^

For instance, to use these scripts on a TORQUE based cluster, one might update
~/.mvmany.template to the following ::

    #!/bin/bash
    #PBS -N $jobname
    #PBS -l nodes=1
    #PBS -l ppn=1
    #PBS -l mem=$memory
    #PBS -l walltime=$walltime
    #PBS -e $logpath/$jobname.e
    #PBS -o $respath/$jobname.txt

    cd $pwd

    $body

Please note that not all SLURM settings have a direct mapping to PBS settings
and that it is up to the user to understand how to properly configure their
cluster job headers.

In general, the user should ensure that each of the variables are properly
defined so that the corresponding values will be written to the final job
scripts. The following variables are replaced based on the job that is being
performed and the parameters passed to the program by the user (or their
default values):

.. tabularcolumns:: |p{6cm}|p{9cm}|
=================================  =============================================
  **Variable**                      **Purpose**
=================================  =============================================
  $jobname                          Unique name for the current job
  $memory (2G)                      Amount of memory to provide each job.
  $walltime (3:00:00)               Define amount of time to be assigned to jobs
  $logpath                          Directory specified for writing logs
  $respath                          Directory sepcified for writing results
  $pwd                              current working dir when mvmany is run
  $body                             Statements of execution
=================================  =============================================

Command Line Arguments
++++++++++++++++++++++

mvmany.py exposes the following additional arguments for use when running
the script.

.. option:: --mv-path PATH

   Set path to mvtest.py if it's not in PATH

.. option:: --logpath PATH

   Path to location of job's error output

.. option:: --res-path PATH
   Path to location of job's results

.. option:: --script-path PATH

   Path for writing script files

.. option:: --template FILENAME

   Specify a template other than the default

.. option:: --snps-per-job INTEGER

   Specify the number of SNPs to be run at one time

.. option:: --mem STRING

   Specify the amount of memory to be requested for each job

.. option:: --wall-time

   Specify amount of time to be requested for each job


The option, --mem, is dependent on the type of input that is being used as well
as configurable options to be used. The user should perform basic test runs
to determine proper settings for their jobs. By default, 2G is used, which is
generally more than adequate for binary pedigrees, IMPUTE and transposed
pedigrees. Others will vary greatly based on the size of the dataset and the
settings being used.

The option, --wall-time, is largely machine dependent but will vary based on
the actual dataset's size and completeness of the data. Users should perform
spot tests to determine reasonable values. By default, the requested wall-time
is 3 days, which is sufficient for a GWAS dataset, but probably not
sufficient for an entire whole exome dataset and the time required will depend
on just how many SNPs are being analyzed by any given node.

In general, mvmany.py accepts all arguments that mvtest.py accepts, with the
exception of those that are more appropriately defined by mvmany.py itself.
These include the following arguments ::

    --chr
    --snps
    --from-bp
    --to-bp
    --from-kb
    --to-kb
    --from-mb
    --to-mb

To see a comprehensive list of the arguments that mvmany.py can use simply
ask the program itself ::

    mvmany.py --help



Users can have mvmany split certain types of jobs up into pieces and can
specify how many independent commands to be run per job. At this time,
mvmany.py assumes that imputation data is already split into fragments and
doesn't support running parts of a single file on multiple nodes.

The results generated can be manually merged once all nodes have completed
execution.