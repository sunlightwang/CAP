# CAPTRE
This is the script repository for CAPTRE analysis (cap profiling for translational efficiency).

The repository includes a few sections. 

1. Identification of TSS-isoform pairs with significant divergent TE (bootstraping + FDR control) 
   * altTSS_sigDivTE
2. List of 5p-UTRs isforms without alternative splicing 
   * 3T3_5UTR_CS
3. Sequence features responsible for TE regulation 
   * Seq_feat

##Installation: 
* CAPTRE requires
  - [Perl] (https://www.perl.org/) >=5.10 with
    - [Math-Random-OO] (http://search.cpan.org/~dagolden/Math-Random-OO-0.22/) >=0.22
  - [R] (https://www.r-project.org/) >=3.2

* Download scripts

  `git clone https://github.com/sunlightwang/CAPTRE.git`

  and `cd` to each directory to execute 

##Contact:
Xi Wang (xi dot wang at mdc-berlin dot de)

##Citation:
[1] Xi Wang\*, Jingyi Hou\*, Claudia Quedenau, Wei Chen (2016) Pervasive isoform-specific translational regulation via alternative transcription start sites in mammals. Under revision. 
