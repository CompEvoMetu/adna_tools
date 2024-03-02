# adna_tools

============================================================

INDEX:

 1) How to set the paths for Python scripts, aDNA Tools, and ped-sim
 2) How to install the aDNA Tools Python wheel package and ped-sim
 3) How to filter or interpolate VCF files given a genetic map
 4) How to generate random individuals
 5) How to select specific samples from the 1000G dataset
 6) How to compute the minor allele frequency values
 7) How to split a VCF file by extracting a specific chromosome
 8) How to filter a VCF file according to a minor allele frequency threshold
 9) How to create a script to simulate multiple instances of pedigrees
10) How to compute (and plot if necessary) the dynamic allele sharing coefficient and the non-normalized mismatch coefficient among samples in a VCF file
11) Split the VCF file into multiple files (one for each sample)
12) How to pseudo-haploidize a VCF file
13) How to normalize the mismatch coefficient
14) How to append multiple VCF files
15) How to insert the ngsRelate tag with the frequency information in a VCF file
16) How to split a VCF file with missing values into multiple VCF files for subsequent dynamic ASC analysis
17) How to check a VCF file for compatibility with aDNA Tools
18) How to merge multiple individuals from multiple input files into a newly created VCF file
19) How to estimate MAF and alternate frequency values from multiple VCF files

============================================================

1) HOW TO SET THE PATHS FOR PYTHON SCRIPTS, ADNA TOOLS, AND PED-SIM:

   Define the path for aDNA Tools and for the in-house scripts (use your specific folders):

   $ ADNA_PATH="/path/to/python/packages/adna_tools/"; export ADNA_PATH                                                                  
   $ IN_HOUSE_PYTHON_PATH="/path/to/python/packages/"; export IN_HOUSE_PYTHON_PATH                                                                          
   $ PEDSIM_PATH="/path/to/PedSim"; export PEDSIM_PATH

   These lines can be added to the .profile (or .bashrc) file so that at every log-in the variables are automatically created (e.g., cd ~; nano .profile).
   Note that you have to install your own version of PedSim and specify its path to the PEDSIM_PATH variable (see last paragraph in 2). Also, you have to create the    
   /python/packages/ directory in your specific path.

============================================================

2) HOW TO INSTALL THE ADNA TOOLS PYTHON WHEEL PACKAGE AND PEDSIM:

   ADNA TOOLS

   $ git clone https://github.com/CompEvoMetu/adna_tools.git
   
   After cloning the repository, location of aDNA Tools package should be: /path/to/python/packages/adnaTools-2021.3.22-py3-none-any.whl
   
   $ pip3 install --target=$IN_HOUSE_PYTHON_PATH adnaTools-2022.5.14-py3-none-any.whl --upgrade

   The main available functions are:
   • Filter/Interpolate VCF files given genetic maps [3]
   • Generate random individuals given minor allele frequencies [4]
   • Compute minor allele frequencies given a VCF file [6]
   • Split a VCF file by extracting a chromosome [7]
   • Filter a VCF file given an MAF threshold [8]
   • Create a script to simulate multiple instances of pedigrees [9]
   • Create plots on how the allele-sharing coefficient varies among samples in a VCF file

   PEDSIM
   Note: Be sure you are on a mounted disk which allows you to execute "executable" files.

   $ git clone https://github.com/williamslab/ped-sim.git                                                                                                                              
   $ cd ped-sim                                                                                                                                                                        
   $ make

============================================================

3) HOW TO FILTER OR INTERPOLATE VCF FILES GIVEN A GENETIC MAP:

   Examples:

   • Documentation:
     $ python $ADNA_PATH/filter_vcf.py --help

   • Filter VCF and genetic map with both autosomes and X chromosome:
     $ python $ADNA_PATH/filter_vcf.py --vcf tuscany.vcf --prefix interpolated_tuscanians
                                       --f_maps ./genetic_map/female_chr{}.txt --m_maps ./genetic_map/male_chr{}.txt --x_chr 23

   • Filter VCF and genetic map with only autosomes:
     $ python $ADNA_PATH/filter_vcf.py --vcf tuscany_autosomes.vcf --prefix interpolated_tuscanians
                                       --f_maps ./genetic_map/female_chr{}.txt --m_maps ./genetic_map/male_chr{}.txt

   • Filter VCF and genetic map with only X chromosome:
     $ python $ADNA_PATH/filter_vcf.py --vcf tuscany_x.vcf --prefix interpolated_tuscanians
                                       --f_maps ./genetic_map/female_chr{}.txt --m_maps ./genetic_map/male_chr{}.txt --chrs [] --x_chr 23

   • Interpolate genetic map values given VCF with both autosomes and X chromosome:
     $ python $ADNA_PATH/filter_vcf.py --vcf tuscany.vcf --prefix interpolated_tuscanians
                                       --f_maps ./genetic_map/female_chr{}.txt --m_maps ./genetic_map/male_chr{}.txt --x_chr 23 --interpolate

   • Interpolate genetic map values given VCF with only autosomes:
     $ python $ADNA_PATH/filter_vcf.py --vcf tuscany_autosomes.vcf --prefix interpolated_tuscanians
                                       --f_maps ./genetic_map/female_chr{}.txt --m_maps ./genetic_map/male_chr{}.txt --interpolate

   • Interpolate genetic map values given VCF with only X chromosome:
     $ python $ADNA_PATH/filter_vcf.py --vcf tuscany_x.vcf --prefix interpolated_tuscanians
                                       --f_maps ./genetic_map/female_chr{}.txt --m_maps ./genetic_map/male_chr{}.txt --chrs [] --x_chr 23 --interpolate

============================================================

4) HOW TO GENERATE RANDOM INDIVIDUALS:

   Output(s), varying in accordance with the given parameters:
   • VCF file with the randomized individuals
   • ASC file with the allele sharing coefficient estimates between all individuals
   • LOG file with the feedback on the randomization process
   • FAM file with the sexes of the generated individuals

   Note: The frequency file must contain the alternate allele frequency values and not the minor allele frequency (MAF) values.
         It is erroneously named as MAF but the alternate allele frequency values must be passed.

   Examples:

   • Documentation:
     $ python $ADNA_PATH/rnd_pop.py --help

   • Generate random individuals with both autosomes and X chromosome and sexes given a MAF file:
     $ python $ADNA_PATH/rnd_pop.py --ref maf_template.frq --prefix random_template --n_male 3 --n_female 2 --x_chr X

     Note: The parameter x_chr could be omitted as the default value identifying X chromosome SNPs is 'X'.

   • Generate random individuals with autosomes (ignoring X chromosome as no sex specifications are given) with an MAF file as 
     input:
     $ python $ADNA_PATH/rnd_pop.py --ref maf_template.frq --prefix random_template --n 10

   • Generate random individuals with both autosomes and X chromosome and sexes given a VCF file with MAF specified with a tag 
     (within the VCF file):
     $ python $ADNA_PATH/rnd_pop.py --ref tuscany.vcf --maf_ref tag --prefix random_template --n_male 3 --n_female 2 --x_chr 23

   • Generate random individuals with autosomes (ignoring X chromosome as no sex specifications are given) with a VCF file and 
     MAF defined via tag (within the VCF file):
     $ python $ADNA_PATH/rnd_pop.py --ref tuscany.vcf --maf_ref tag --prefix random_template --n 10 --x_chr 23

     Note: The parameter x_chr is specified as the VCF file contains also X chromosome data when found since the sexes are 
           not required, it will be ignored.
           If the VCF file does not contain X chromosome data then there is no need to specify this parameter.

   • Generate random individuals with autosomes and sexes with a VCF file and MAF defined via tag (within the VCF file):
     $ python $ADNA_PATH/rnd_pop.py --ref tuscany_autosomes.vcf --maf_ref tag --prefix random_template --n_female 5 --n_male 9

   • Generate random individuals with autosomes, with no sex specified, given a VCF file and MAF defined via tag (within the 
     VCF file):
     $ python $ADNA_PATH/rnd_pop.py --ref tuscany_autosomes.vcf --maf_ref tag --prefix random_template --n 10

   • Generate random individuals with X chromosome and sexes given a VCF file and MAF defined via tag (within the VCF file):
     $ python $ADNA_PATH/rnd_pop.py --ref tuscany_x.vcf --maf_ref tag --prefix random_template --n_female 5 --n_male 9 --x_chr 23

   • Generate random individuals with both autosomes and X chromosome and sexes given MAF file and specifying column indices:
     $ screen python $ADNA_PATH/rnd_pop.py --ref TSI_maf_01.frq --prefix random_TSI --n_male 1 --n_female 1
                                           --x_chr X --chr_idx 0 --pos_idx 1 --id_idx 2 --ref_idx 3 --alt_idx 4 --maf_idx 5

============================================================

5) HOW TO SELECT SPECIFIC SAMPLES FROM THE 1000G DATASET

   The filter_id.sh takes a  text file as input—containing  IDs of the desired samples—and loops across
   all the 1000G files, initially to filter for bi-allelic SNPs and, later on, to retain only the  IDs
   that were specified in the input text file:

   $ screen -m bash ./filter_id.sh TSI_ID.txt

   Input:
   [TSI_ID.txt] contains the list of the IDs of interest.

   Output:
   [TSI_bi_snps.vcf] will contain bi-allelic SNPs of only Tuscany individuals.
   Additional files are created during the script, they can be removed afterward if not useful
   (more details on those files below).

   IMPORTANT REMINDER:

   The filtered files are stored in the directory where the script is launched.
   If different parameters are needed, you should update the script accordingly.
   1000G files directory should be set in the filter_id.sh as well as the bcftools path manually.

   • source: variable for the 1000G dataset location.
   • cmd: variable for the bcftools command.
   • output: variable for the output prefix file(s) (following bi-allelic SNP filtering).
   • list_txt: file that will contain the list of the newly filtered VCF files which will
     be merged (automatically created by the script).
   • merge_vcf: variable identifying the VCF file created following the merging.
   • filtered_vcf: variable identifying the final VCF file containing only bi-allelic SNPs
     for Tuscany individuals.

============================================================
 
6) HOW TO COMPUTE THE MINOR ALLELE FREQUENCY VALUES
   
   Given an  input VCF file  ('vcf_file'),  vcf2maf computes minor  allele frequencies for all SNPs
   and stores them in a given output file ('maf_file'). Unrecognized SNPs samples are stored in a
   log file ('vcf_file'.log). Only entries '0' (+1  reference base), '1' (+1 alternate base), '0|0'
   (+2 reference base),  '0|1'  (+1 reference base,  +1 alternate base)   '1|0' (+1 reference base,
   +1 alternate base) and '1|1'  (+2  alternate base)  are considered.  If all SNPs in the VCF file
   have the aforementioned format, the log file will be empty.  In the given example are  shown
   phased samples but also unphased files are accepted and automatically recognized (however it can
   be explicitly declared via the parameter phased, see also the documentation).  If missing values
   are encountered, they are simply ignored and the count for the frequency is based on the sum of
   the found reference and alternate values.
   
   Example:

   • Documentation:
     $ python3 $ADNA_PATH/vcf2maf.py --help

   • $ screen python3 $ADNA_PATH/vcf2maf.py --vcf_file TSI_bi_snps.vcf --maf_file TSI_maf.frq

============================================================

7) HOW TO SPLIT A VCF FILE BY EXTRACTING A SPECIFIC CHROMOSOME

   Given a VCF file and a chromosome ID as inputs, it creates two new VCF files, where one contains only the specified
   chromosome and one with all the remaining chromosomes. By default, the X chromosome ('X') is extracted.  Optionally,
   for the X chromosome,  there is a pseudo-haploid conversion of the format,  e.g., '0' becomes '0|0' and '1' becomes
   '1|1'.

   Example:

   • Documentation:
     $ python3 $ADNA_PATH/extract_chr.py --help

   • $ screen python3 $ADNA_PATH/extract_chr.py --vcf_file TSI_bi_snps.vcf --prefix 'TSI_bi_snps' --chr_id 23 --pseudo-haploid True

============================================================

8) HOW TO FILTER A VCF FILE ACCORDING TO A MINOR ALLELE FREQUENCY THRESHOLD

   Given  an input VCF or MAF file and a minor allele frequency threshold value, it filters all the SNPs with MAF
   smaller  than  the given threshold.  The frequency values are assumed  to be stored within the VCF file with a
   specific tag (default 'AFngsrelate') or, if an MAF file is passed, within a specific column. A MAP file can also
   be passed to maintain the consistency of the SNPs between MAF and MAP files.

    Examples (ADNA_PATH="/mnt/.../python/packages/adna_tools/"; export ADNA_PATH):

    • Documentation:                                                                                                                                                                   
      $ python $ADNA_PATH/maf_filter.py --help

    • Filter VCF file given threshold:                                                                                                                                                 
      $ python $ADNA_PATH/maf_filter.py --maf_file pops.vcf --prefix filter_pops --maf_threshold 0.05
   
    • Filter FRQ file given threshold:                                                                                                                                                 
      $ screen python3 $ADNA_PATH/maf_filter.py --maf_file pops_maf.frq --prefix filter_pops --maf_threshold 0.1

============================================================

9) HOW TO CREATE A SCRIPT TO SIMULATE MULTIPLE INSTANCES OF PEDIGREES

   To create a Unix shell script to generate all the desired pedigrees you can use the make command from the aDNA Tools.

   Examples (ADNA_PATH="/mnt/.../python/packages/adna_tools/"; export ADNA_PATH):

    • Documentation:
        $ python $ADNA_PATH/make.py --help

    • Initialize a script to create multiple instances (10) of a set of pedigrees (defined in the pedigrees.cfg):
        $ python $ADNA_PATH/make.py --config ./def/pedigrees.cfg --n 10 --prefix founders --gmap inter_biSNPs_MAF_01_TSI.map
                                    --maf ./MAF/TSI_maf_01.frq --m_founders 10 --f_founders 10
                                    --working_folder ./def --export init_pop.sh

   It needs a  configuration file  ('config')  and the number of samples  ('n') for each relationship specified in the
   configuration file.  Minor allele frequency ('maf') and genetic map ('gmap') files, as well as the desired number of
   male ('m_founders') and female ('f_founders') founders, must also be specified. Optionally, the location of DEF and
   IDX  files  can be specified (by default the files are assumed to be  in the current folder). Other relevant (yet
   optional) parameters include the tag identifying the X chromosome ('x') and the name of the exported  script ('export')
   and log file   ('log_file',  in which all the combinations  of the given  pedigrees/relationships  are summarized).
   Notably, the creation of the script assumes that the variables $PEDSIM_PATH and $ADNA_PATH are defined and identify
   ped_sim and adnaTools folders respectively. By default, the interference map ('int_map') is assumed to be located in
   the ped-sim folder (i.e., $PEDSIM_PATH/interfere/nu_p_campbell_X.tsv). The configuration file contains a list of the
   relationships of interest identified by their IDX files.
   Ex.: 'pedigrees.cfg':
   #par-off_2nd-cous_inb.idx
   #par-off_1st-cous_inb.idx
   #par-off_sibs_inb.idx
   #siblings.idx
   #parent-offspring.idx
   grandparent-grandchild.idx
   grandparent-grandchild_sibs_inb.idx
   grandparent-grandchild_1st-cous_inb.idx
   grandparent-grandchild_2nd-cous_inb.idx
   avuncular.idx
   avuncular_1st-cous_inb.idx
   avuncular_2nd-cous_inb.idx
   half-siblings.idx
   half-siblings_sibs_inb.idx

   Notably, the lines starting with # are considered as comments and therefore ignored.

   • The abstract form of a DEF file:
     A DEF file is a definition file used by ped-sim to create pedigrees. It contains a compact form defining a
     pedigree of interest (see ped-sim website for more details on its syntax). An abstract DEF file is a
     non-complete DEF file that will be used by this script to generate all multiple combinations of sexes within a
     specified pedigree.
     Ex.: Parent-Offspring pedigree - 'parent-offspring.def'

     def parent-offspring 1 2
     1 1 2 {}
     2 1 1 1:1_2 {}

     Whereas the general syntax is analogous to the one of ped-sim, there are two sets of brackets '{}' that act as
     placeholders and will be used by this function to format different combinations of this pedigree (according
     to specific sexes; this syntax mimics the 'print('text {} text').format(variable)' adopted by Python).

     In this specific example, we have two possible combinations for this pedigree when we account for the sexes of
     the individuals. Specifically, for the first combination, generation 1 will have one male and one female
     individuals (i.e., the parents) and a female individual in generation 2 (i.e., the offspring of the individuals
     in generation 1). Alternatively, the second combination, generation 1 will still have one male and one female
     individuals (i.e., the parents) and a male individual in generation 2 (i.e., the offspring of the individuals
     in generation 1).

     Thus, this script will generate two combinations of the defined pedigree considering the different sexes the
     each individual may assume.
     The details on the possible sexes that can be assigned to each individual are provided by the IDX file.

     Additional examples of abstract DEF and IDX files are available at '/mnt/NAS/projects/2019_igor/TSI/def'.

   • The IDX file:
     The IDX file should specify all the possible sexes each individual can have. This script will be responsible then
     to generating all possible combinations. Beyond that, it should also provide details about all the relationships of
     interest within the specific pedigree.

     Ex.: Parent-Offspring pedigree - 'parent-offspring.idx'
     NOTE: The DEF and IDX files MUST have the same name.

     <M=1 F=0>
     <G=1:[(0,1)]>
     <G=2:[0 1]>
     <IDX=0:(1,1)>
     <IDX=1:(1,2)>
     <IDX=2:(2,1)>
     #IDX1	IDX2 RELATIONSHIP DEGREE
     0	2	Parent-Offspring	1
     1	2	Parent-Offspring	1

     Notably:
     '<M=1 F=0>' simply clarifies the association between used IDs and sexes.

     '<G=1:[(0,1)]>' refers to generation 1 in the DEF file (i.e., G=1) and provides the list of acceptable sexes
     for the individuals in generation one, i.e., in this case, the two parents, only one accepted combination (0,1).
     Notably, (1,0) would only be redundant and thus can be ignored.

     '<G=2:[0 1]>' refers to the generation 2 in the DEF file (i.e., G=2). Once again, it provides the list of accepted
     sexes for the single individual, i.e., female 0 and male 1.

     Notably, if the DEF file has a generation with two or more branches/individuals the possible combinations are
     enclosed in parentheses: (0,0). A generation with a single branch/individual would not need parentheses: 0.
     To identify all possible combinations, the accepted sexes are enclosed by squared brackets:
     [(0,0) (0,1) (1,0) (1,1)] or [0 1].

     In certain cases, it may be of interest to ignore a generation (i.e., its sex combinations) in the context of the
     symmetrical evaluation. This can be indicated by adding an asterisk (*) next to the generation of interest in the
     IDX file as follows:
     '<G=1*:[(0,1)]>' i.e., the generation one won't influence the analysis of symmetrical combination.

     The IDX lines identify the indices of the individuals that will be stored in the ped-sim generated VCF file:
     '<IDX=0:(1,1)>' the individual stored in column 0 corresponds to the one in generation 1 and branch 1 of the
     DEF file.
     '<IDX=1:(1,2)>' the individual stored in column 1 corresponds to the one in generation 1 and branch 2 of the
     DEF file.
     '<IDX=2:(2,1)>' the individual stored in column 2 corresponds to the one in generation 2 and branch 1 of the
     DEF file.

     Notably, its general form is <IDX=n:(generation,branch)> where 'n' is the n-th column identifying the individual
     in the ped-sim generated VCF file (with 0 being the first individual) and the tuple (generation,branch)
     identifying the generation and branch of said individual in the DEF file.
     These associations are necessary for the script to be able to retrieve the sexes of each individual of interest
     from all possible combinations.

     The last part of the file is structured as a tabbed table:
     #IDX1	IDX2 RELATIONSHIP DEGREE
     0	2	Parent-Offspring	1
     1	2	Parent-Offspring	1

     '#IDX1	IDX2 RELATIONSHIP DEGREE': the first line representing its header.
     '0	2	Parent-Offspring	1' and
     '1	2	Parent-Offspring	1': The first and second columns identify the individuals to be considered for the
     relationship specified in the third column. The fourth column defines the degree of said relationship.
     Thus, in this case, we have two parent-offspring relationships between individuals 0 and 2 (first entry) as well as
     1 and 2 (second entry).
     Notably, the genders of those relationships vary for each combination and this script will list them in the
     logging file.

     Additional examples of abstract DEF and IDX files are available at '/mnt/NAS/projects/2019_igor/TSI/def'.

   • The logging file:
       - Four columns whose names are specified in the header:
           1) RELATIONSHIP: Name of the relationship of interest inclusive of the sexes of the individuals.
              Ex.: Half-siblings_F_F (two females half-siblings)
                   Half-siblings_M_F (male-female half-siblings)
                   Half-siblings_M_M (two male half-siblings)
          2) DEGREE: The degree of the relationships.
              Ex.: 2
          3) (COMBINATION,IDX1,IDX2): A list of triplets that allow  to identify the  individuals  (possessing the
             specified relationship) within the generated VCF files. Specifically, the first integer identifies the
             combination number whereas the second and third indices identify the columns  within the VCF file of said
             individuals.

             Ex.: [(1, 0, 1), (1, 2, 3), (4, 0, 1), (4, 2, 3)], i.e.,

             (1, 0, 1): for the combination 1 (whose file name is determined by the REFERENCE column - e.g.,
             half-siblings_{}.vcf - and the combination number, thus, in the first case 'half-siblings_1.vcf'), the
             individuals identified by indices 0 (the first sample column of the VCF file) and 1 (the second sample
             column of the VCF file) are half-siblings.

             (1, 2, 3): for the combination 1 (whose file name is determined by the REFERENCE column - e.g.,
             half-siblings_{}.vcf - and the combination number, thus, in this case, 'half-siblings_1.vcf'), the
             individuals identified by indices 2 and 3 are half-siblings.

             (4, 0, 1): for the combination 4 (whose file name is determined by the REFERENCE column - e.g.,
             half-siblings_{}.vcf - and the combination number, thus, in this case, 'half-siblings_4.vcf'), the
             individuals identified by the indices 0 and 1 are half-siblings.

             (4, 2, 3): for the combination 4 (whose file name is determined by the REFERENCE column - e.g.,
             half-siblings_{}.vcf - and the combination number, thus, in this case, 'half-siblings_4.vcf'), the
             individuals identified by indices 2 and 3 are half-siblings.

     4) REFERENCE: The general syntax to identify the files associated with each mentioned relationship.
        Ex.: half-siblings_{}.vcf, i.e., The individuals generated with the combination 1 are stored in
        half-siblings_1.vcf. In a more general way, the individuals generated with the combination N are stored in
        half-siblings_N.vcf (where N is an integer).

============================================================

10) HOW TO COMPUTE (AND PLOT IF NECESSARY) THE DYNAMIC ALLELE SHARING COEFFICIENT AND THE NON-NORMALIZED MISMATCH COEFFICIENT AMONG SAMPLES IN A VCF FILE

    Given a VCF file, it computes the dynamic allele sharing coefficient and the non-normalized mismatch coefficient (given 
    specific window length and shift) among all samples in the file for a given set of chromosomes (one figure and logging
    file per chromosome will be generated). The window is applied over the SNPs list.
    If you want to print the indices instead of the IDs stored in the VCF file (or other specified labels) you can use the
    parameter --use_labels set to False. The flag 'plots' allows to creation of output figures if necessary (--plots True).
    The coefficients are computed and the plots are created together with the log files where all the used values are stored.

    NOTE 1: The non-normalized mismatch coefficient is still not finalized and it just computes the ratio between the number
    of observed mismatches and the total number of checks.

    NOTE 2: Among the options you can specify the chromosomes of interest. There is however a caveat. If your VCF file contains
    all chromosomes, by only selecting a few, the function will still analyze all the files (and hence all chromosomes) and
    only later will filter the results retaining solely the one of interest. Thus, whenever you are interested in only a few
    chromosomes, the best approach would be (in order to speed up the process) to filter the VCF file and retain only the
    chromosomes you are interested in.

    NOTE 3: The parameter 'phased' is an optional parameter that by default is set to None and lets the script determine
    if the available data in the VCF are phased or unphased (the automatic decision is made by simply looking at structure of
    the first available SNP: if the separator '|' is found the data are assumed phased; if the separator '/' is found the data
    are assumed unphased). However, it is possible to explicitly specify the type of data by setting the 'phased' parameter to
    True (phased) or False (unphased).

    Example:

    • Documentation:
      $ python3 $ADNA_PATH/dynamic_kinship.py --help

    • Compute dynamic ASC on all chromosomes (each chromosome has one figure and one logging file) with a non-overlapping 
      window and no plots:
      $ screen python3 $ADNA_PATH/dynamic_kinship.py --vcf_file samples.vcf --prefix asc_samples
                                                     --win_length 500000 --win_shift 500000 --plots False

    • Compute dynamic ASC on only the X chromosome with a non-overlapping window:
      $ screen python3 $ADNA_PATH/dynamic_kinshio.py --vcf_file samples.vcf --prefix asc_samples
                                                     --win_length 500000 --win_shift 500000 --chromosomes X

    • Compute dynamic ASC on a series of chromosomes (each chromosome has one figure and one logging file) with a non- 
      overlapping window:
      $ screen python3 $ADNA_PATH/dynamic_kinship.py --vcf_file samples.vcf --prefix asc_samples
                                                     --win_length 500000 --win_shift 500000 --chromosomes [1,2,3,X]

    • Compute dynamic ASC on a series of chromosomes (each chromosome has one figure and one logging file) with 50% 
      overlapping window:
      $ screen python3 $ADNA_PATH/dynamic_kinship.py --vcf_file samples.vcf --prefix asc_samples
                                                     --win_length 500000 --win_shift 250000 --chromosomes [1,2,3,X]

    • Compute dynamic ASC on a series of chromosomes (each chromosome has one figure and one logging file) with 50% 
      overlapping window
      with only selected relationships and associated labels:
      $ screen python3 $ADNA_PATH/dynamic_kinship.py --vcf_file samples.vcf --prefix asc_samples
                                                     --win_length 500000 --win_shift 250000 --chromosomes [1,2,3,X]
                                                     --idx_rel [0,1,2] --labels ["Relation0","Relation1","Relation2"]

============================================================

11) SPLIT VCF FILE IN MULTIPLE FILES (ONE FOR EACH SAMPLE)

    Given an input VCF file with M samples (or a text file with a list of VCF files, one per line) it creates M VCF
    files each with a single sample. The newly created files will have the name of the source followed by an index.

    Example:

    • Documentation:
      $ python3 $ADNA_PATH/split_vcf.py --help

    • Split a VCF file:
      $ python3 $ADNA_PATH/split_vcf.py --source population.vcf

      Output files for the population.vcf with M samples:
      population1.vcf
      population2.vcf
      ...
      populationM.vcf

    • Split a list of VCF files whose names are stored in a text file:
      $ python3 $ADNA_PATH/split_vcf.py --source list.txt

============================================================

12) HOW TO PSEUDO-HAPLOIDIZE A VCF FILE:
    
    To apply the pseudo-haploidization to a VCF file you can call the haploidize.py function.
    The function will parse the VCF input file and create one or multiple copies (each with unique randomized haploidizations)
    stored with the same file name with an additional suffix (by default '_haploid') and index (to identify the different v 
    haploidization copies).
    Whereas the only required parameter is the input VCF file, there are other parameters that can be specified (see the 
    documentation).

    Example:

    • Documentation:
      $ python3 $ADNA_PATH/haploidize.py --help

    • Single run of haploidization for a VCF file (the output VCF file will contain 0s and 1s):
      $ python3 $ADNA_PATH/haploidize.py --vcf_file sample.vcf

    • Three runs of haploidization for a VCF file (the output VCF files will contain single 0 and 1 values):
      $ python3 $ADNA_PATH/haploidize.py --vcf_file sample.vcf --n 3

    • Two runs of haploidization for a VCF file (the output VCF files will contain the two-element values 0|0 or 1|1):
      $ python3 $ADNA_PATH/haploidize.py --vcf_file sample.vcf --n 2 --sigle False

============================================================

13) HOW TO NORMALIZE THE MISMATCH COEFFICIENT:

    When applying the dynamic_asc function (see also paragraph 10) over a VCF file of interest, not only the allele sharing
    coefficient will be computed but also the mismatch coefficient. The latter is non-normalized, i.e., the estimated value
    does not account for the baseline mismatch coefficient pertaining to unrelated samples. Thus, following its computation
    there are additional steps that need to be performed in order to obtain a normalized set of estimates.
    The first step would be to compute the dynamic mismatch coefficient for a VCF file containing two or more unrelated
    individuals. This can be simply done by executing the dynamic_asc function (with the same parameters used to compute the
    coefficients on the samples of interest). As an example, we can consider a VCF file containing a set of unrelated founders:
  
    • Example:  
      $ python3 $ADNA_PATH/dynamic_asc.py --vcf_file founders_1.vcf --prefix unrelated --win_length 5000 --win_shift 5000

    The mismatch coefficient log file will then contain the estimates (computed among all unrelated samples) for each of the
    analyzed windows.
    
    • See generated output:
      $ less -S unrelated_1_MSM.log

    In this case, the index 1 indicates the output for chromosome 1.

    Upon these values, we are going to compute the baseline mismatch coefficient for each window by calling the function
    msm_baseline:

    • Documentation:
      $ python3 $ADNA_PATH/msm_baseline.py --help

    • Example:
      $ python3 $ADNA_PATH/msm_baseline.py --log_file unrelated_1_MSM.log --prefix msm_baseline_median --mode median

    The function considers a LOG file pertaining to the mismatch coefficient (see dynamic_asc documentation for more details)
    and computes the related baseline (according to the indicated mode) by considering all sample combinations stored in the
    LOG file. Thus, if you apply the dynamic_asc function over a VCF file containing unrelated individuals, the produced LOG
    file will have the window-averaged mismatch coefficient for all combinations of samples. By passing that file to the
    msm_baseline one can estimate the baseline reference for each window. The mode parameter indicates how to compute this 
    baseline.
    By default, the mode="median" will consider the median value within each window. Alternative implemented modes are "mean", 
    which computes the average mismatch value, "max" the highest observed estimate, or "all" where multiple columns, one for 
    each mode, are exported.
    The output file will be a text file (CVS formatted) with the extension BSL:

    • See generated output:
      $ less -S msm_baseline_median.bsl

    Whereas the first and second columns will contain the window IDs and the sample size respectively, the third column will 
    store the baseline value as requested (with its column named accordingly, e.g., median). If the chosen mode was "all",   
    then columns three, four, and five would store respectively the median, mean, and maximum baseline estimates.

    Once the baseline file(s) is (are) computed, with it (them) you can normalize all the relevant mismatch coefficient LOG 
    files by calling the normalize_msm function:

    • Documentation:
      $ python3 $ADNA_PATH/normalize_msm.py --help

    • Example:
      $ python3 $ADNA_PATH/normalize_msm.py --log_file samples_1_MSM.log --bsl_file msm_baseline_median.bsl --prefix 
      normalized_samples_1

    The function will normalize the values stored in the LOG file and will create a new file with the given prefix and an 
    automatic suffix identifying the applied normalization:

    • See generated output:  
      $ less -S normalized_samples_1_BSL_median.log

============================================================

14) HOW TO APPEND MULTIPLE VCF FILES:

    To append multiple VCF files (e.g., if you have one VCF file per chromosome) you can call the function append_vcf:
    append_vcf(txt_file: str, prefix: str):
    Given a text file containing a list of VCF files, it generates a new VCF file by appending each
    VCF following the order defined in the input file (one filename per line). No check on the
    consistency of the different VCF files is done, lines are read and written to the output file
    unchecked.

    NOTE: Only the commented lines of the first VCF file are retained.

    :param str txt_file: Text file containing the list of VCF files to be merged (appended) one after the
    other.
    :param str prefix: Output filename
    Basically, you create a text file where each line contains a VCF filename (ex.: list.txt) then you pass that file, 
    together with a prefix (ex.: merged_ancient) for the output VCF file that will be generated, to the function.
    
    • Documentation:
      $ python3 $ADNA_PATH/append_vcf.py --help

    • Example:
      $ python3 $ADNA_PATH/append_vcf.py --txt_file filenames_list.txt --prefix merged_filename

============================================================

15) HOW TO INSERT THE NGSRELATE TAG WITH THE FREQUENCY INFORMATION IN A VCF FILE:

    To insert the ngsRelate tag with the frequency information you can use the command add_AF_tag,
    you need to pass the input VCF file and the frequency file (one frequency value per line and no
    other information). The assumption is that the first line contains the frequency value for the
    first SNP in the VCF file (and so on).

    • Documentation:
      $ python3 $ADNA_PATH/add_AF_tag.py --help

    • Example:
      $ python3 $ADNA_PATH/add_AF_tag.py --vcf_file ancient.vcf --af_file ancient_frequency.frq

============================================================

16) HOW TO SPLIT A VCF FILE WITH MISSING VALUES INTO MULTIPLE VCF FILES FOR SUBSEQUENT DYNAMIC ASC ANALYSIS:

    When creating multiple VCF files from a given input VCF file with missing values a new folder
    is created and in each new file the available common SNPs of a couple of individuals are stored.
    Thus, given an input VCF file (ancient.vcf) with 4 individuals (0, 1, 2, 3), 6 new files
    will be created (inside a new folder "ancient"): ancient_0_1.vcf (shared SNPs between
    individual 0 and 1), ancient_0_2.vcf (shared SNPs between individual 0 and 2),
    ancient_0_3.vcf, ancient_1_2.vcf, ancient_1_3.vcf, ancient_2_3.vcf.

    • Documentation:
      $ python3 $ADNA_PATH/ancient_vcf.py --help

    • Example:
      $ python3 $ADNA_PATH/ancient_vcf.py --vcf_file ancient.vcf

============================================================

17) HOW TO CHECK A VCF FILE FOR COMPATIBILITY WITH ADNA TOOLS

    The script fix_vcf, allows to check the compatibility of a given VCF file with our in-house scripts.
    It controls chromosome/position consistency as well as supported SNPs values (0s and 1s) while also
    filtering out missing values (if necessary). Beyond filtering it can also fix chromosome/position
    problems. Basically, you can pass a single VCF file or a text file containing a list of VCF files
    (one filename per line) to be examined. It evaluates each input for compatibility issues and if necessary
    it generates new VCF files (with the additional given prefix) where all misplaced entries are rearranged
    and duplicate entries as well as unsupported SNP lines are discarded.
    
    • Documentation:
      $ python3 $ADNA_PATH/fix_vcf.py --help

    • Check multiple VCF files passed as a text list while also retaining SNPs with at most 2 missing values per sample:
      $ python3 $ADNA_PATH/fix_vcf.py --filename list.txt --missing_values 2 --prefix fixed_ --feedback False

    • Check a single VCF file while also filtering all SNPs with at least one missing value (i.e. retain SNPs with 0 missing 
      values):
      $ python3 $ADNA_PATH/fix_vcf.py --filename tuscany.vcf --missing_values 0 --prefix fixed_ --feedback False

============================================================

18) HOW TO MERGE MULTIPLE INDIVIDUALS FROM MULTIPLE INPUT FILES INTO A NEWLY CREATED VCF FILE

    The script merge_vcf allows to creation of new VCF output from a given text file containing a list of VCF files.
    The created file contains the individuals from each file. Positions that are not found in certain samples
    are marked with missing values. Common SNPs with conflicting reference and/or alternate values are automatically
    discarded (if necessary a log file provides all the details).

    • Documentation:
      $ python3 $ADNA_PATH/merge_vcf.py --help

    • Example:
      $ python3 $ADNA_PATH/merge_vcf.py --txt_file "list.txt" --prefix "merged_data"

============================================================

19) HOW TO ESTIMATE MAF AND ALTERNATE FREQUENCY VALUES FROM MULTIPLE VCF FILES

    The script vcf2freq allows to estimation of MAF and alternate frequency values from a single VCF file or from a text file 
    where each line identifies a VCF file to consider. Alongside these values, the script also generates LOG and ERR files 
    with the number of available samples at each SNP position for each  observed value (LOG listing all the details, ERR only 
    details pertaining to discarded SNPs).
    Four output files will be created: {prefix}MAF.frq, {prefix}ALT.frq, {prefix}FRQ.log, {prefix}FRQ.err for MAF values, 
    alternate frequency values, log, and error files respectively. Each FRQ file will contain three columns: CHR (chromosome), 
    POS (SNP position), FRQ (frequency value) unless the flag 'only_frq' is set to True, in which case only frequency values 
    will be saved (by default 'only_frq' is set to False).
    The LOG and ERR files will contain eight columns: CHR (chromosome), POS (SNP position), REF (the reference allele), ALT 
    (the alternate allele),
    N_REF (total number of samples for the reference allele), N_ALT (total number of samples for the alternate allele), STATUS 
    (indicating if the SNP was retained in the final output, i.e., OK (only LOG), or it was discarded due to 
    reference/alternate/non-bi-allele inconsistencies or unknown values, INFO (additional information, used to provide the 
    list of observed SNPs and their occurrences in discarded samples).
    When a SNP is discarded the REF, ALT, N_REF, N_ALT columns are filled with '-'.

    • Documentation:
      $ python3 $ADNA_PATH/vcf2freq.py --help

    • Example:
      $ python3 $ADNA_PATH/vcf2freq.py --txt_file "list.txt" --prefix "neo_data"

============================================================
