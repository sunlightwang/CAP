# CAPTRE
This is the script repository for CAPTRE analysis (cap profiling for translational efficiency).

The repository includes scripts for the following tasks. 

1. Identification of TSS-isoform pairs with significant divergent TE (bootstraping + FDR control) 

2. Construction of 5p-UTRs isforms without alternative splicing 

3. Sequence features responsible for TE regulation 

##Installation: 
* CAPTRE requires
  - [Perl] (https://www.perl.org/) >=5.10 with
    - [Scalar-List-Utils] (http://search.cpan.org/~pevans/Scalar-List-Utils-1.45/) >=1.45
    - [Math-Random-OO] (http://search.cpan.org/~dagolden/Math-Random-OO-0.22/) >=0.22
  - [R] (https://www.r-project.org/) >=3.2 with 
    - [Hmisc] (https://cran.r-project.org/web/packages/Hmisc/index.html)
  - [bedtools] (http://bedtools.readthedocs.io/en/latest/) >=2.17.0

* Download scripts

  `git clone https://github.com/sunlightwang/CAPTRE.git`

  then `cd` to each directory and refer to README.sh inside to run. Three folders for the three tasks are listed below:
  * altTSS_divTE
  * 3T3_5UTR_CS
  * Seq_feat

##Contact:
Xi Wang (xi dot wang at mdc-berlin dot de)

##Citation:
[1] Xi Wang\*, Jingyi Hou\*, Claudia Quedenau, Wei Chen (2016) Pervasive isoform-specific translational regulation via alternative transcription start sites in mammals (http://msb.embopress.org/content/12/7/875). Molecular Systems Biology, 12: 875.

