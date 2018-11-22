SimulateCNVs

Maintainer: Yue "July" Xing
Author: Yue "July" Xing
Version: 1.0
Date: 06/27/2018


Description:
A tool for simulating CNVs for WES or WGS data. It simulates rearranged genomes, short reads (fastq) and bam files in a single command. There are several ways and distributions to choose from to generate desired CNVs.


Installation:
No installation required. Simply type the following and it will be ready to use: git clone https://github.com/YJulyXing/SimulateCNVs.git
Or manually download the source codes at https://github.com/YJulyXing/SimulateCNVs.

Requirements:
1. General use: Python 2.7 (https://www.python.org/download/releases/2.7/). Required python packages: argparse, random, os, subprocess, math, sys, time
2. To generate short reads (fastq) outputs (see requirements for ART_illumina (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)):
  GNU g++ 4.0 or above (http://gcc.gnu.org/install)
  GNU gsl library (http://www.gnu.org/s/gsl/)
3. To generate bam outputs:
  Samtools (http://samtools.sourceforge.net/)
  BWA (http://bio-bwa.sourceforge.net/)
  picard 2.15.0 (https://broadinstitute.github.io/picard/)
  GATK (https://software.broadinstitute.org/gatk/)


Usage:
SimulateCNVs.py [-h] -Type {g,e} -G GENOME_FILE [-T TARGET_REGION_FILE]
                [-e_cnv TARGET_CNV_LIST] [-e_chr TARGET_CNV_CHR] [-e_tol TARGET_CNV_TOL] [-e_cl TARGET_CNV_LEN_FILE]
                [-o_cnv OUT_CNV_LIST] [-o_chr OUT_CNV_CHR] [-o_tol OUT_CNV_TOL] [-o_cl OUT_CNV_LEN_FILE]
                [-ol OVERLAP_BPS] [-g_cnv GENOME_CNV_LIST] [-g_chr GENOME_CNV_CHR] [-g_tol GENOME_CNV_TOL]
                [-g_cl GENOME_CNV_LEN_FILE] [-em] [-min_len CNV_MIN_LENGTH] [-max_len CNV_MAX_LENGTH]
                [-min_cn MIN_COPY_NUMBER] [-max_cn MAX_COPY_NUMBER] [-p PROPORTION_INS] [-f MIN_FLANKING_LEN]
                [-ms {random,uniform,gauss}] [-ml {random,uniform,gauss,user}] [-c COVERAGE]
                [-fs FRAG_SIZE] [-s STDEV] [-l READ_LENGTH] [-tf TARGET_REGION_FLANK] [-pr]
                [-q_min MIN_BASE_QUALITY] [-q_max MAX_BASE_QUALITY] [-clr CONNECT_LEN_BETWEEN_REGIONS]
                [-o OUTPUT_DIR] [-rn REARRANGED_OUTPUT_NAME] [-n NUM_SAMPLES] [-sc]
                [-ssr] [-sb] [-picard PATH_TO_PICARD] [-GATK PATH_TO_GATK]


Arguments:
1. Optional arguments:

|   Parameter               | Default |    Explanation                             | Restrictions                    |
                               value
| :-----------------------: | :-----: | :----------------------------------------- | :------------------------------ |
| -h, --help                | -       | show this help message and exit            | -                               |

2. Mandatory arguments:

|   Parameter               | Default |    Explanation                             | Restrictions                    |
                               value
| :-----------------------: | :-----: | :----------------------------------------- | :------------------------------ |
| -Type {g,e}               | -       | simulation for WGS or WES                  | -                               |
| -G GENOME_FILE            | -       | Reference genome FASTA file                | -                               |

3. Arguments for simulating rearranged genomes for WES data:

|   Parameter               | Default |    Explanation                             | Restrictions                    |
                               value
| :-----------------------: | :-----: | :----------------------------------------- | :------------------------------ |
| -T TARGET_REGION_FILE     | -       | Target region file                         | Must be used and can only be    |
                                                                                      used with WES simulation.      
| -e_cnv TARGET_CNV_LIST    | -       | A user-defined list of CNVs overlapping    | One and only one of -e_cnv,     |
                                        with target regions                           -e_chr, -e_tol and -e_cl can   
                                                                                      be used with WES simulation    
                                                                                      to generate CNVs overlapping   
                                                                                      with target regions. If -e_cnv is       
                                                                                      provided, -em, -f, -ms, -ml,   
                                                                                      -ol, -min_cn, -max_cn,         
                                                                                      -min_len and -max_len will be  
                                                                                      ignored for CNVs overlapping   
                                                                                      with target regions.                    
| -e_chr TARGET_CNV_CHR       | -       | Number of CNVs overlapping with target     | Same as above.                  |
                                          regions to be generated on each chromosome            
| -e_tol TARGET_CNV_TOL       | -       | Total number of CNVs overlapping with      | Same as above.                  |
                                          target regions to be generated across the  
                                          genome (an estimate)
| -e_cl TARGET_CNV_LEN_FILE   | -       | User supplied file of CNV length for CNVs  | Must be used with -ml user.     |
                                          overlapping with target regions                Can’t be used with -o_cnv, 
                                                                                         -o_chr and -o_tol. Otherwise 
                                                                                         same as above.
| -o_cnv OUT_CNV_LIST       | -         | A user-defined list of CNVs outside of     | One and only one of -o_cnv,     |
                                         target regions                                 -o_chr, -o_tol and -o_cl can 
                                                                                        be used with WES simulation   
                                                                                        to generate CNVs outside of 
                                                                                        target regions. If -o_cnv is provided,
                                                                                        -em, -f, -ms, -ml, -ol, 
                                                                                        -min_cn, -max_cn, -min_len and -max_len will be ignored for 
                                                                                        CNVs outside of target regions.
| -o_chr OUT_CNV_CHR        | -       | Number of CNVs outside of target regions to be      | Same as above.           |
                                         generated on each chromosome
| -o_tol OUT_CNV_TOL        | -       | Total number of CNVs outside of target regions to   | Same as above.           |
                                         be generated across the genome (an estimate)  
| -o_cl OUT_CNV_LEN_FILE    | -       | User supplied file of CNV length for CNVs           | Must be used with -ml user.|
                                         outside of target regions                             Can’t be used with -e_cnv, 
                                                                                               -e_chr and -e_tol.  
                                                                                               Otherwise same as above.
| -ol OVERLAP_BPS           | 100     | For each CNV overlapping with target regions,       | Can only be used with WES  |
                                         number of minimum overlapping bps                     simulation.

4. Arguments for simulating rearranged genomes for WGS data:

|   Parameter               | Default |    Explanation                             | Restrictions                    |
                               value  
| :-----------------------: | :-----: | :----------------------------------------- | :------------------------------ |
| -g_cnv GENOME_CNV_LIST    | -       | A user-defined list of CNVs outside of     | One and only one of -g_cnv, 
                                        target regions                                -g_chr, -g_tol and -g_cl can 
                                                                                      be used with WGS simulation to generate CNVs.
| -g_chr GENOME_CNV_CHR     | -       | Number of CNVs overlapping with target     | Same as above.                  |
                                        regions to be generated on each chromosome 
| -g_tol GENOME_CNV_TOL     | -       | Total number of CNVs overlapping with      |  Same as above.                 |
                                        target regions to be generated across the 
                                        genome (an estimate) 
| -g_cl GENOME_CNV_LEN_FILE | -       | User supplied file of CNV length           | Same as above.                  |

5. General arguments for simulating rearranged genomes with CNVs:

|   Parameter               | Default |    Explanation                             | Restrictions                    |
                               value  
| :-----------------------: | :-----: | :----------------------------------------- | :------------------------------ |
| -em                       | -       | Exclude missing sequences for CNV          | -                               |
                                         simulation
| -min_len CNV_MIN_LENGTH   | 1000    | Minimum CNV length in bps                  | -                               |
| -max_len CNV_MAX_LENGTH   | 100000  | Maximum CNV length in bps                  |
| -min_cn MIN_COPY_NUMBER   | 2       | Minimum copy number for insertions         | -                               |
| -max_cn MAX_COPY_NUMBER   | 10      | Maximum copy number for insertions         | -                               |
| -p PROPORTION_INS         | 0.5     | Proportion of insertions                   | -                               |
| -f MIN_FLANKING_LEN       | 50      |  Minimum length between each CNV           | -                               |
| -ms {random,uniform,gauss}| random  | Distribution of CNVs                       | -                               |
| -ml {random,uniform,      | random  | Distribution of CNV length                 | -ml user must be used with      |
      gauss,user}                                                                     -e_cl and/or -o_cl. If -ml user 
                                                                                      is used, -min_len and -max_len 
                                                                                      will be ignored.

6. Arguments for simulating short reads (fastq):

|   Parameter               | Default |    Explanation                             | Restrictions                    |
                               value  
| :-----------------------: | :-----: | :----------------------------------------- | :------------------------------ |
| -c COVERAGE               | 20      | Fold coverage on target regions to be      | -                               |
                                          generated for each genome
| -fs FRAG_SIZE             | 100     | Mean fragment size to be generated         | -                               |
| -s STDEV                  | 20      | Standard deviation of fragment sizes       | -                               |
| -l READ_LENGTH            | 50      | Read length of each short read             | -                               |
| -tf TARGET_REGION_FLANK   | 0       | Length of flanking region up and down      | -tf takes place after -clr      |
                                          stream of target regions to be sequenced    (target regions are first
                                                                                      connected per request of the 
                                                                                      user, and then flanking regions 
                                                                                      up and down stream of target 
                                                                                      regions are included for sequencing). Only works with 
                                                                                      WES simulation.
| -pr                       | -       | Select if paired-end sequencing            | For paired-end sequencing, must 
                                                                                      use the option -pr. If using single-end sequencing, mean 
                                                                                      fragment size (-fs) and 
                                                                                      standard deviation of fragment 
                                                                                      size (-s) will be ignored.
| -q_min MIN_BASE_QUALITY   | 0       | Minimum base quality for short reads       |-                                |
                                         simulation
| -q_max MAX_BASE_QUALITY   | 80      | Maximum base quality for short reads       | -                               |
                                         simulation 
| -clr                      | -       | Maximum length bwtween target regions to   | -tf takes place after -clr      |
 CONNECT_LEN_BETWEEN_REGIONS              connect the target regions                  (target regions are first 
                                                                                      connected per request of the 
                                                                                      user, and then flanking regions 
                                                                                      up and down stream of target 
                                                                                      regions are included for sequencing). Only works with 
                                                                                      WES simulation.

7. Arguments for other simulation parameters:

|   Parameter                | Default value     |    Explanation                     | Restrictions                    | 
| :------------------------: | :---------------: | :--------------------------------- | :------------------------------ |
| -o OUTPUT_DIR              | simulation_output | Output directory                   | Will be generated if not exist. |
| -rn REARRANGED_OUTPUT_NAME | test              | Prefix of the rearranged outputs   | Do not include directory name.  |
| -n NUM_SAMPLES             | 1                 | Number of test samples to be       | Must be >= 1.                   |
                                                    generated 
| -sc                        | -                 | Simulation for control genome      | -                               |
| -ssr                       | -                 | Simulate short reads (fastq) files | If the final output is bam      |
                                                                                         file(s), must first simulate 
                                                                                         short reads (-ssr) and then 
                                                                                         bam (-sb).
| -sb                        | -                 | Simulate bam files                 | Same as above.                  |
| -picard PATH_TO_PICARD     | -                 | Absolute path to picard            | If the final output is bam      |
                                                                                         file(s) (using -ssr and -sb), 
                                                                                         must provide absolute paths to 
                                                                                         picard and GATK.
| -GATK PATH_TO_GATK         | -                 | Absolute path to GATK              | Same as above.                  |


Inputs:
1. Sequence of a reference genome in fasta format. Anything after "_" or " " at header line will be ignored.
2. For WES simulation, a tab delimited file of target regions in the order of chromosome, start and end. Header should not be included. Chromosome name should strictly match chromosome name in (1). Target regions should be sorted in ascending order.
3. For -e_cnv, -o_cnv and g_cnv, a tab delimited file generated by SimulateCNVs, in the order of chromosome, start, end, length and copy number. Header should be included. Chromosome name should strictly match chromosome name in (1).
4. For -e_cl, -o_cl and -g_cl, a tab delimited file of desired CNV lengths in the order of chromosome, CNV length and number of CNV of that length. Header should not be included. Chromosome name should strictly match chromosome name in (1).


Outputs:
1. Rearranged genome(s) (fasta)
Target regions for rearranged genome(s) (bed)
Control genome (fasta, if -sc is chosen)
Target regions for control (bed, always generated in case -clr is chosen to make changes to target regions)
2. List(s) of CNVs overlapping with target regions (bed)
List(s) of CNVs outside of target regions (bed, if chosen to generate CNVs outside of target regions)
3. Short reads for rearranged genome(s) (fastq, if -ssr is chosen)
Short reads for control genome(s) (fastq, if -sc and -ssr is chosen)
4. Indexes for the control genome (dict, fai, sa, etc., if -ssr and -sb is chosen and no indexes exist in the output directory)
5. Bam file(s) and index(es) for rearranged genome(s) (bam and bai, if -ssr and -sb is chosen)
Bam file(s) and index(es) for control genome (bam and bai, if -sc, -ssr and -sb is chosen)


Examples:
Simulation for WGS and WES data is similar. WES simulation, however, need more parameter settings.
a. WGS data simulation:
  1. Simulate 10 CNVs of length 1 kb to 10 kb on each chromosome randomly, and at least 100 bps between each 2 CNVs. Don’t generate CNVs on missing sequences. Make a pair of test and control genomes.

  SimulateCNVs/SimulateCNVs.py -Type g -G <input_fasta> -o <output_dir> \
                                -g_chr 10 -sc -min_len 1000 -max_len 10000 -f 100 -em


  2. Simulate approximately 100 CNVs of length 1 kb to 10 kb on the whole genome, 30% of which are insertions, and at least 200 bps between each 2 CNVs. CNV start points and CNV lengths both follow gauss distribution. Generate CNVs on missing sequences. Make 10 test samples with prefix “test_wgs” and does not make any control. Make bam files as final output. Single-end sequencing is used.

  SimulateCNVs/SimulateCNVs.py -Type g -G <input_fasta> -o <output_dir>  -g_tol 100 -p 0.3 \
                                -min_len 1000 -max_len 10000 -f 200 -ms gauss -ml gauss -ssr -sb \
                                -n 10 -rn test_wgs -picard <absolute_path_to_picard> \
                                -GATK < absolute_path_to_GATK>

  3. Distribution of CNV lengths are user provided. CNV start points follow uniform distribution. The copy numbers range from 5 to 15 for insertions. Don’t generate CNVs on missing sequences. Make a pair of test and control. Make short reads (fastq) file as final output using paired-end sequencing, 40 fold coverage, 100 bp read length, mean fragment size 300 bp and standard deviation of mean fragment size 10.
  SimulateCNVs/SimulateCNVs.py -Type g -G <input_fasta> -o <output_dir> \
                                -ml user -g_cl <CNV length file> -min_cn 5 -max_cn 15 \
                                -em -ms uniform -sc -pr -ssr -c 40 -fs 300 -s 10 -l 100

b. WES data simulation:
  4. Simulate 10 CNVs overlapping with target regions, and 1 CNV outside of target regions randomly on each chromosome using default lengths, copy numbers, minimum distance between each of the 2 CNVs and proportion of insertions. For each CNV overlapping with target regions, the overlapping length is not less than 90 bps. CNV start points and lengths follow gauss distribution. Don’t generate CNVs on missing sequences. Make 5 test samples and control. Generate short reads (fastq) files by default settings, using paired-end sequencing.
  
  SimulateCNVs/SimulateCNVs.py -Type e -G <input_fasta> -T <target_region> -o <output_dir> \
                                -e_chr 10 -o_chr 1 -ol 90 -ms gauss -ml gauss -em -n 5 -sc -pr -ssr

  5. Simulate CNVs overlapping with target regions from the provided CNV list. Simulate approximately 20 CNVs outside of target regions randomly on the whole genome with default settings. For CNVs outside of target regions, don’t generate CNVs on missing sequences. Make a pair of test and control genome.

  SimulateCNVs/SimulateCNVs.py -Type e -G <input_fasta> -T <target_region> -o <output_dir> \
                                -e_cnv <list_of_CNV_overlapping_with_target regions> -o_tol 20 -em -sc 

  6. Simulate approximately 20 CNVs overlapping with target regions on the whole genome, and at least 100 bps between each 2 CNVs. Don’t generate CNVs outside of target regions. Don’t generate CNVs on missing sequences. Paired-end sequencing, with minimum base quality is 20 and maximum base quality is 60. Make a pair of test and control. The final outputs are bam files.
  SimulateCNVs/SimulateCNVs.py -Type e -G <input_fasta> -T <target_region> -o <output_dir> \
                                -e_tol 20 -f 100 -em -sc -pr -q_min 20 -q_max 60 -ssr -sb \
                                -picard <absolute_path_to_picard> -GATK <absolute_path_to_GATK>

  7. Simulate CNVs overlapping with target regions and outside of target regions from provided files of CNV lengths. If the length between 2 target regions are smaller than 100 bps, connect them as 1 target region. Don’t generate CNVs on missing sequences. Make 10 test samples and control. Use paired-end sequencing; sequence 50 bp up and down stream of the target regions (after connecting the target regions) as well. The final output is short reads (fastq) files with coverage of 40. 

  SimulateCNVs/SimulateCNVs.py -Type e -G <input_fasta> -T <target_region> -o <output_dir> \
                                -ml user -e_cl <length_file_1> -o_cl <length_file_2> \
                                -clr 100 -em -n 10 -sc -pr -tf 50 -f 40 -ssr 


=====================================================================================================================

ReplaceNs.py

Maintainer: Yue "July" Xing
Author: Yue "July" Xing
Version: 1.0
Date: 06/27/2018


Description:
A small program to fix genomes which have too many ‘N’s to generate desired CNVs. It replaces all ‘N’s in the genome sequence to ‘A’s, ‘T’s, ‘G’s, or ‘C’s randomly.


Installation:
It is included in the package of SimulateCNVs.


Requirements:
Python 2.7 (https://www.python.org/download/releases/2.7/)


Usage:
ReplaceNs.py [-h] -i INPUT_FASTA_FILE -o OUTPUT_FASTA_FILE


Arguments:
-h: help
-i: input genome sequence in fasta format
-o output genome sequence in fasta format
