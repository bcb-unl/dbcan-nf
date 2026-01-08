nf-core/dbcanmicrobiome: Usage
================================

.. warning::
   Please read this documentation on the nf-core website: https://nf-co.re/dbcanmicrobiome/usage

   Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files.

Introduction
------------

This pipeline supports three analysis modes for CAZyme annotation in microbiome data:

- **Short reads** (``--type shortreads``): Assembly-based analysis for Illumina short-read data using MEGAHIT
- **Long reads** (``--type longreads``): Assembly-based analysis for PacBio/Nanopore long-read data using Flye
- **Assembly free** (``--type assemfree``): Direct annotation without assembly using seqtk and DIAMOND blastx

The assembly free mode is particularly useful for large datasets where assembly is computationally expensive or when you want to avoid potential assembly artifacts. It directly converts reads to FASTA format and uses DIAMOND blastx to search against the CAZyme database, followed by abundance calculation using dbcan_asmfree.

Samplesheet input
-----------------

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

.. code-block:: bash

   --input '[path to samplesheet file]'

Multiple runs of the same sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``sample`` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

.. code-block:: csv
   :name: samplesheet.csv

   sample,fastq_1,fastq_2
   CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
   CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
   CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz

Full samplesheet
~~~~~~~~~~~~~~~~

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where ``TREATMENT_REP3`` has been sequenced twice.

.. code-block:: csv
   :name: samplesheet.csv

   sample,fastq_1,fastq_2
   CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
   CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
   CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
   TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
   TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
   TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
   TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,

+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Column   | Description                                                                                                                                                                    |
+==========+======================================================================================================================================================================================+
| ``sample`` | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (``_``). |
+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``fastq_1`` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``fastq_2`` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
+----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

An `example samplesheet <../assets/samplesheet.csv>`_ has been provided with the pipeline.

Running the pipeline
--------------------

The typical command for running the pipeline is as follows:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome --input ./samplesheet.csv --outdir ./results  -profile docker

This will launch the pipeline with the ``docker`` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

.. code-block:: bash

   work                # Directory containing the nextflow working files
   <OUTDIR>            # Finished results in specified location (defined with --outdir)
   .nextflow_log       # Log file from Nextflow
   # Other nextflow hidden files, eg. history of pipeline runs and old logs.

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a ``yaml`` or ``json`` file via ``-params-file <file>``.

.. warning::
   Do not use ``-c <file>`` to specify parameters as this will result in errors. Custom config files specified with ``-c`` must only be used for `tuning process resource specifications <https://nf-co.re/docs/usage/configuration#tuning-workflow-resources>`_, other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome -profile docker -params-file params.yaml

with:

.. code-block:: yaml
   :name: params.yaml

   input: './samplesheet.csv'
   outdir: './results/'
   <...>

You can also generate such ``YAML``/``JSON`` files via `nf-core/launch <https://nf-co.re/launch>`_.

Updating the pipeline
~~~~~~~~~~~~~~~~~~~~~

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

.. code-block:: bash

   nextflow pull nf-core/dbcanmicrobiome

Reproducibility
~~~~~~~~~~~~~~~

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the `nf-core/dbcanmicrobiome releases page <https://github.com/nf-core/dbcanmicrobiome/releases>`_ and find the latest pipeline version - numeric only (eg. ``1.3.1``). Then specify this when running the pipeline with ``-r`` (one hyphen) - eg. ``-r 1.3.1``. Of course, you can switch to another version by changing the number after the ``-r`` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse `parameter files <#running-the-pipeline>`_ to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

.. tip::
   If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

Core Nextflow arguments
-----------------------

.. note::
   These options are part of Nextflow and use a *single* hyphen (pipeline parameters use a double-hyphen)

``-profile``
~~~~~~~~~~~

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

.. important::
   We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from `https://github.com/nf-core/configs <https://github.com/nf-core/configs>`_ when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the `nf-core/configs documentation <https://github.com/nf-core/configs#documentation>`_.

Note that multiple profiles can be loaded, for example: ``-profile test,docker`` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If ``-profile`` is not specified, the pipeline will run locally and expect all software to be installed and available on the ``PATH``. This is *not* recommended, since it can lead to different results on different machines dependent on the computer environment.

- ``test``
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- ``docker``
  - A generic configuration profile to be used with `Docker <https://docker.com/>`_
- ``singularity``
  - A generic configuration profile to be used with `Singularity <https://sylabs.io/docs/>`_
- ``podman``
  - A generic configuration profile to be used with `Podman <https://podman.io/>`_
- ``shifter``
  - A generic configuration profile to be used with `Shifter <https://nersc.gitlab.io/development/shifter/how-to-use/>`_
- ``charliecloud``
  - A generic configuration profile to be used with `Charliecloud <https://hpc.github.io/charliecloud/>`_
- ``apptainer``
  - A generic configuration profile to be used with `Apptainer <https://apptainer.org/>`_
- ``wave``
  - A generic configuration profile to enable `Wave <https://seqera.io/wave/>`_ containers. Use together with one of the above (requires Nextflow ``24.03.0-edge`` or later).
- ``conda``
  - A generic configuration profile to be used with `Conda <https://conda.io/docs/>`_. Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

Short reads modes: subsample & coassembly
------------------------------------------

Only relevant when ``--type shortreads``.

Subsample mode
~~~~~~~~~~~~~~

- **Purpose**: downsample each sample before assembly to reduce compute and to quickly sanity-check the pipeline.
- **Parameters**:
  - ``--subsample``: enable subsampling.
  - ``--subsample_size``: number of reads per file to keep (default ``20000000`` in config).
- **Behavior**:
  - Applies per-sample before MEGAHIT, using ``seqtk sample``.
  - Mutually exclusive with ``--coassembly``.
- **Example**:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --type shortreads \
     --input samplesheet.csv \
     --outdir results_subsample \
     --subsample \
     --subsample_size 5000000 \
     -profile docker

Coassembly mode
~~~~~~~~~~~~~~~

- **Purpose**: co-assemble all short-read samples together to improve contig continuity and shared feature detection.
- **Parameters**:
  - ``--coassembly``: enable coassembly across all samples.
- **Requirements & behavior**:
  - Needs at least 2 samples; pipeline will error if fewer.
  - Combines all reads (paired/single preserved) and runs a single MEGAHIT assembly.
  - Mutually exclusive with ``--subsample``.
- **Example**:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --type shortreads \
     --input samplesheet.csv \
     --outdir results_coassembly \
     --coassembly \
     -profile docker

``-resume``
~~~~~~~~~~

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see `this blog post <https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html>`_.

You can also supply a run name to resume a specific run: ``-resume [run-name]``. Use the ``nextflow log`` command to show previous run names.

``-c``
~~~~~~

Specify the path to a specific config file (this is a core Nextflow command). See the `nf-core website documentation <https://nf-co.re/usage/configuration>`_ for more information.

Custom configuration
--------------------

Resource requests
~~~~~~~~~~~~~~~~~

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified `here <https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18>`_ it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the `max resources <https://nf-co.re/docs/usage/configuration#max-resources>`_ and `tuning workflow resources <https://nf-co.re/docs/usage/configuration#tuning-workflow-resources>`_ section of the nf-core website.

Custom Containers
~~~~~~~~~~~~~~~~~

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the `biocontainers <https://biocontainers.pro/>`_ or `bioconda <https://bioconda.github.io/>`_ projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the `updating tool versions <https://nf-co.re/docs/usage/configuration#updating-tool-versions>`_ section of the nf-core website.

Custom Tool Arguments
~~~~~~~~~~~~~~~~~~~~~

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the `customising tool arguments <https://nf-co.re/docs/usage/configuration#customising-tool-arguments>`_ section of the nf-core website.

nf-core/configs
~~~~~~~~~~~~~~~

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the ``nf-core/configs`` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the ``-c`` parameter. You can then create a pull request to the ``nf-core/configs`` repository with the addition of your config file, associated documentation file (see examples in `nf-core/configs/docs <https://github.com/nf-core/configs/tree/master/docs>`_), and amending `nfcore_custom.config <https://github.com/nf-core/configs/blob/master/nfcore_custom.config>`_ to include your custom profile.

See the main `Nextflow documentation <https://www.nextflow.io/docs/latest/config.html>`_ for more information about creating your own configuration files.

If you have any questions or issues please send us a message on `Slack <https://nf-co.re/join/slack>`_ on the `#configs channel <https://nfcore.slack.com/channels/configs>`_.

Running in the background
--------------------------

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow ``-bg`` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use ``screen`` / ``tmux`` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

Nextflow memory requirements
----------------------------

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in ``~/.bashrc`` or ``~./bash_profile``):

.. code-block:: bash

   NXF_OPTS='-Xms1g -Xmx4g'

