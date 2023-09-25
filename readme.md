
# DTy (dcgMLST typer)
nomenclature engine for distributed core genome MLST scheme (dcgMLST) and pre-defined HierCC clustering


# INSTALLATION:

DTy was developed and tested in Python >=3.8. It depends on several Python libraries: 
~~~~~~~~~~
click
numba
numpy
pandas
biopython
pyarrow
fastparquet
~~~~~~~~~~

All libraries can be installed using pip: 

~~~~~~~~~~
pip install click numba numpy pandas biopython pyarrow fastparquet
~~~~~~~~~~
DTy also calls two 3rd party programs:

~~~~~~~~~~
ncbi-blast+
diamond
~~~~~~~~~~

Both can be installed via 'apt' in UBUNTU:
~~~~~~~~~~
sudo apt install -y ncbi-blast+ diamond
~~~~~~~~~~

The whole environment can also be installed in conda:


~~~~~~~~~~
conda create --name dty python==3.11
conda activate dty
conda install -c conda-forge biopython numba numpy pandas click pyarrow fastparquet
conda install -c bio-conda blast diamond
~~~~~~~~~~

The installation process normally finishes in <10 minutes. 


# Quick Start (with examples)
## Get allelic and HierCC callings
~~~~~~~~~~~
$ cd /path/to/DTy/
$ python DTy.py -q examples/GCF_005221305.fna -o examples/GCF_005221305.alleles
~~~~~~~~~~~
## Generate tab-delimited allelic profiles
~~~~~~~~~~~
python alleles2profile.py -p examples/examples examples/*.alleles
~~~~~~~~~~~

The whole calculation finishes in <2 minutes with 8 CPU threads. 


# USAGE:
## DTy.py - allelic callings and HierCC clusters & species predictions

~~~~~~~~~~~~~~
$ python DTy.py --help
Usage: DTy.py [OPTIONS]

Options:
  -q, --query TEXT        input assembly. fasta or fastq format, can be
                          gziped.
  -o, --outfile TEXT      alleles output file
  -n, --n_thread INTEGER  n_threads [default:8]
  --help                  Show this message and exit.
~~~~~~~~~~~~~~~~~

## alleles2profile.py - summarization of multiple allelic files into one profile file. 

~~~~~~~~~~~~~
$ python alleles2profile.py --help
Usage: alleles2profile.py [OPTIONS] [FNS]...

Options:
  -c, --scheme TEXT  scheme name. default: Neisseria_dcgMLST
  -p, --prefix TEXT  prefix for parquet/tsv output of the allelic profiles
                     [required]
  --help             Show this message and exit.
~~~~~~~~~~~~~~~~

# Outputs:
## DTy.py
### DTy.py generates screen outpus:

~~~~~~~~~~~~~
$ python DTy.py -q examples/GCF_005221305.fna -o examples/GCF_005221305.alleles
examples/GCF_005221305.fna      HC1130.1050.760.400.300.100.50.20.10.5.2=159.5768.5768.5768.5768.5768.N.N.N.N.N species=Neisseria subflava      reference=GCA_005221305 allelic_distance=80
~~~~~~~~~~~~~

that gives out HC & species predictions. It also gives out evidence for the predictions (closest reference and distance to it). 

### DTy.py also generates an allelic calling results (examples/GCF_005221305.alleles here) like:

~~~~~~~~~~~~~
$ head -4 examples/GCF_005221305.alleles
>CCD84_RS00010 value_md5=f7d49f149ad872340b86448e3602f565 CIGAR=CCD84_RS00010_1:264M39D48M9I12M12I96M3D18M9I1089M accepted=2 reference=GCF_005221305.fna identity=0.846 coordinates=NZ_CP039887.1:1..1569
ATGACGTTAGCTGAATTTTGGCCGCTGTGCCTTCGCCGCCTTCACGAAATGTTGCCTGCCGGGCAGTTTGCGCAATGGATTGCGCCTTTGACCGTGGGCGAAGAAAACGGCGTATGGGTGGTGTATGGTAAAAACCAATTTGCCTGCAATATGCTCAAAAGCCAGTTTGCCGCCAAAATTGACGCCGTGCGTGCCGAATTAGTGCCTCAGCAGGCTGCTTTTGCGTTTAAGCCGGGCGTAGGTACGCATTATGAAATGGCGGCTCAGACTGTTGCGCCGGTGCAAGTGCAAGAGGTCATTGAAGTTGAAGAGTGTGTAGAGCCTGTTCAAATGCCTTTGCAAACTGCTGCGCCAATGGAAGAAAATAGGCCGTCTGAAACGGTTTCCAAACCTGCAGCTGCCATGACGGCTGCCGAGATTTTGGCGCAACGCATGAAAAACCTGCCGCATGAGCCTCAAGTGCAAACTACTGCTTCGGCTGAATCTAAAGCAGTTGCCAAAGCCAAAACCGACGCGCAACACGATGCGGAAGAAGCGCGCTACGAACAGACCAATCTGTCGCGTGACTATACATTTGAAACTTTGGTGGAAGGTAAGGGCAACCGCCTTGCTGCCGCAGCCGCCCAGGCGATTGCTGAAAATCCGGGGCAGGGCTACAACCCGTTTTTCTTATACGGCAGTACCGGTTTGGGTAAAACCCACTTGGTGCAAGCCATCGGCAACGAATTGCTGAAAAACCGTCCTGATGCCAAAGTGCGCTATATGCACTCGGATGACTATATCCGCAGCTTTATGAAGGCGGTGCGCAACAATACTTACGATGTATTCAAGCAACAATACAAACAATATGACCTGCTGATTATCGACGATATCCAGTTCATCAAAGGCAAAGACCGTACGATGGAAGAATTCTTTTATCTGTACAACCATTTTCACAATGAGAAAAAACAACTTATCCTGACGTGCGACGTATTGCCTGCCAAAATCGAAGGTATGGATGACCGTCTCAAATCCCGTTTCTCATGGGGTTTGACTTTGGAACTCGAACCGCCCGAATTGGAAATGCGCGTGGCGATTTTGCAGAAAAAGGCAGAAGCGGCCGGTATCAGTATCGAAGACGAAGCCGCTCTGTTTATCGCCAATCTGATCCGTTCCAATGTGCGTGAGCTGGAAGGCGCGTTCAACCGCGTCAGCGCCAGCAGCCGTTTTATGAACCGTCCTGTCATTGACATGGATTTGGCGCGTACGGCTTTGCAGGATATTATTGCCGAAAAACACAAAGTCATTACCGCCGACATCATCATCGATGCGACAGCCAAATACTACCGTATTAAAATCAGTGATATATTGGGCAAAAAACGTACGCGCAATATTGCCCGTCCGCGCCAAGTTGCCATGAGCTTGACCAAAGAGCTGACCATGCTCAGCCTGCCTTCTATCGGTGATGCCTTTGGCGGTCGCGATCACACGACTGTGATGCACGGTGTCAAAGCGGTGGCGAAACTGCGCGAAGAAGATCCCGAATTGGCGCAAGACTACGAAAAACTGCTGATTTTGATTCAGAACTGA
>CCD84_RS00015 value_md5=e5e0896a24a85d02324abe1d0e1085e7 CIGAR=CCD84_RS00015_1:1104M accepted=2 reference=GCF_005221305.fna identity=0.894 coordinates=NZ_CP039887.1:1676..2779
ATGCTGATTTTACAAGCCGATCGCGACAGTCTGCTCAAGCCGTTGCAAGCCGTTACCGGTATTGTCGAACGTCGCCATACTCTGCCGATTTTGTCTAATGTGTTGCTGGAAAGCAAAGACGGACAAACCAAACTTTTGGCAACCGACTTGGAAATCCAAATCAATACCGCCGGCCCTGAAAGTCAGGCAGGCGATTTCCGTATTACGACCAACGCTAAAAAATTCCAAGACATCCTGCGTGCTCTGCCTGACAGTGCGCTGGTGTCACTGGATTGGGCGGACAACCGTTTGACTCTGCGCGCGGGAAAATCCCGTTTTGCCCTGCAAACCTTGCCGGCTGAAGACTTTCCGTTGATGAGCGTCGGCAGCGACGTCAGTGCGACTTTCTCACTGACTCAAGAAACCTTCAAAACCATGCTTTCGCAAGTGCAATACAGCATGGCAGTTCAAGATATCCGCTATTACCTCAACGGTTTGCTGATGCAGGTTGAAGGTAATCAGCTGCGCCTTGTTGCAACCGACGGCCACCGCCTTGCTTATGCGGCCAGTCAAATTGAAGCAGAACTGCCGAAAACGGAAGTGATCCTGCCGCGTAAAACGGTATTGGAACTCTTCAAGCTGTTGAATAATCCGTCCGAGTCCATCACCGTTGAGCTTTTGGACAATCAAGTACGCTTCCAATGCAATGGCACAACCATTGTCAGCAAAGTCATCGACGGCAAGTTCCCTGACTTTAACCGCGTGATTCCTTTGGATAATGACAAGATTTTCCTCGTATCCCGTACCCAGCTTTTGGGTGCACTCGAGCGTGCCGCCATTCTTGCCAATGAAAAATTCCGCGGCGCACGACTGTTCCTGCAGCCTGGTTTGCTGAGTGTCGTATGTAGCAACAACGAGCAGGAAGAAGCGCGCGAAGAGCTGGAAATCGCTTACCAAGGCGGAGAACTCGAAGTCGGTTTCAACATCGGCTACCTGATGGATGTGTTGCGTAACATCCACTCCGACGATATGCAGCTTGCTTTCGGCGATGCCAACCGTTCAACGCTGTTTACTATGCCGAACAATCCGAACTTCAAATACATCGTAATGCCGATGCGTATTTAA
~~~~~~~~~~~~~

The header of each allelic sequence includes MD5 values, identity to reference alleles and coordinates in the query assembly. It also includes a CIGAR that describes the alignments. 


## alleles2profile.py
alleles2profile.py generates two files. 
<prefix>.tsv contains a tab-delimited allelic profile (one strain per line). 
<prefix>.parq contains the same information except that in parquet format, with reduced storage space. 


# Citation and Reproduction Instructions
### Citation for the pipeline and the Neisseria dcgMLST scheme
Ling Zhong, Menghan Zhang, Libing Sun, Yu Yang, Bo Wang, Haibing Yang, Qiang Shen, Yu Xia, Jiarui Cui, Hui Hang, Yi Ren, Bo Pang, Xiangyu Deng, Yahui Zhan, Heng Li, and Zhemin Zhou, Distributed genotyping and clustering of Neisseria strains reveals continual emergence of epidemic meningococcus over a century, submitted

### Reproduction Instructions
All data required for reproduction of the analysis were distributed in this repository under
https://github.com/ADSGF203com/DTy/tree/master/db/Neisseria

These includes:
* references.fas.gz - reference alleles for all pan genes (for calling new alleles)
* cgmlst.genes - A list of core genes used in the dcgMLST scheme
* profile.parq - Allelic profiles of all ~70,000 genomes in parquet format, and can be read using the Pandas library (https://pandas.pydata.org/docs/reference/api/pandas.read_parquet.html). 
* HierCC.tsv.gz - A tab-delimited table consisting of HierCC results for all ~70,000 genomes
* hc.species - A mapping table that specifies correlations between HC1130 and HC1050 and the Neisseria species

