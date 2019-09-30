# FunFinder_Pangenome
#### Data analysis of Funannotate and Orthofinder outputs to Pangenome results

## Summary

This script was written to produce many of the results and figures in Wyka et al. (in preparation). It takes output files from [Funannotate](https://github.com/nextgenusfs/funannotate) and the Orthogroups.txt output from [Orthofinder](https://github.com/davidemms/OrthoFinder) (or similar programs such as OrthoMCL and SiLiX) and conducts pangenomic analyses and produces results in the form of text files and figures.

Results will consist of: (Note: these figures are from mock data)

Pangenome stats
![Pangenome stats](https://github.com/PlantDr430/images/blob/master/Pangenome_stats.png)

Pangenome curve - with power law regression fit
![Pangenome curve - with power law regression fit](https://github.com/PlantDr430/images/blob/master/Pangenome_curve.png)

Pangenome fluidity - follows [pangenome_fluidity](https://github.com/PlantDr430/CSU_scripts/blob/master/pangenome_fluidity.py)
![Pangenome fluidity - computed with bootstraps (resampling)](https://github.com/PlantDr430/images/blob/master/Species_9_pangenome_fluidity.png)

Protein length boxplots for pangenome category (lettering for significance with Kruskal-Wallis Test)
![Protein lengths](https://github.com/PlantDr430/images/blob/master/Protein_lengths.png)

Bar charts for all analyses specified at run time showing proportion of clusters in each pangenomic category if >= 50% of species represented in the cluster have at least 1 protein with the respective annotation (i.e. if 5 or more genomes, out of 10 total genomes, have 1 protein in a cluster that has hits for a secreted protein then that cluster will be classifed as a secreted protein cluster) (this can be altered at run time with -p/--percent) (letters for significance with Fischer Exact test)
![Bar charts](https://github.com/PlantDr430/images/blob/master/secretome_pangenome_bar.png)

Gene Ontology Enrichment Analyses for all types (Biological Process, Molecular Function, and Cellular Component) and all Pangenome categories (Core, Accessory, and Singleton)
![GOEA](https://github.com/PlantDr430/images/blob/master/GOEA_Accessory_MF.png)

## Other output files

All figures as well as text files associated with images for reproduction purposes. Within each isolate directory in the species_results directory there will be separate fasta files containing proteins from each analysis provided, as well as, overall_results file for each analysis (Note: These overall_results for each analysis will not always match the pangenome clusters for each analysis as we have a default 50% cutoff for clusters to be classified as a given annotation category). Final output file will be a large tab deliminated file that can be opened in excel which will contain all clusters their proteins and annotations. 

## Dependencies 

1. [Python 3](https://www.python.org/downloads/)   
    * Scipy
    * Numpy
    * Pandas
    * Matplotlib
    * [goatools](https://github.com/tanghaibao/goatools) <- need to pip install
    * [find_enrichment.py](https://github.com/tanghaibao/goatools/blob/master/scripts/find_enrichment.py)   
    
2. [EffectorP_2.0](http://effectorp.csiro.au/software.html)   
3. go-basic.obo   
`wget http://geneontology.org/ontology/go-basic.obo`
4. goslim_generic.obo   
`wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo`

## Usage

```
usage: ./FunFinder_Pangenome.py [options] -d directory -o output

    Data analysis of Funannotate and Orthofinder outputs to Pangenome results.

optional arguments:
  -h, --help                      show this help message and exit
  -d , --directory                Directory containing folders for processing
  -a  [ ...], --analyses  [ ...]  Analyses to compare [default: all]. Can put "none" or individual analyses separated by commas
  -o , --out                      New output folder. Do not use input directory as output directory
  -f_p , --fischer_pvalue         Pvalue cutoff for Fischer exact test (annotation bar charts) [default: 0.05]
  -k_p , --kruskal_pvalue         Pvalue cutoff for Kruskal-Wallis test (protein length boxplot) [default: 0.05]
  -b_p , --benjamini_pvalue       Pvalue cutoff for Benjamini-Hochberg in GO Enrichment Analysis [default: 0.05]
  -al , --alpha                   Test-wise alpha for GO Enrichment Analysis [default: 0.05]
  -p , --percent                  Percent of isolates within a cluster containing proteins with analysis hit, to consider the cluster as such [default: 0.50]
  -c , --cpus                     Number of cores to use for multiprocessing of fluidity [default: 1]
  -max_sub , --max_subsamples     Max number of subsamples to run on N genomes sampled for fluidity. [default: 50000]
  --max_off                       Turn off the max subsamples. This will cause the script sample ALL possible combinationsfor N genomes
  --NO_GO                         Do not perform GO enrichment analysis if iprscan is in analyses [default: ON]
  --EFFECTORP_PATH                Path to EffectorP.py if not set in $PATH
  --ENRICHMENT_PATH               Path to find_enrichment.py if not set in $PATH
  --GOBASIC_PATH                  Path to go-basic.obo if not in current run directory
  --GOSLIM_PATH                   Path to goslim_generic.obo if not in current run directory
```

## Optional analyses:
tmhmm
phobius
signalp
effectors
merops
dbcan
pfam
iprscan
antismash

## Files 

1. Multiple annotations files from annotate_misc Funannotate output folder:

      * annotations.antismash.clusters.txt   
      * annotations.antismash.txt    
      * annotations.dbCAN.txt   
      * annotations.genes-products.txt   
      * annotations.iprscan.txt   
      * annotations.merops.txt   
      * annotations.pfam.txt   
      
          For example the annotations.pfam.txt would look like this:
          
          Spec1_9583-T1  db_xref	PFAM:PF00389   
          Spec1_8343-T1  db_xref	PFAM:PF00004
 
      * phobius.results   
      * signalp.results   
      * tmhmm.results - computed seperately as Funannotate currently doesn't run this   
      
          These will be the direct short output formats of each of these runs (Phobius, SignalP, and TMHMM). I used these over the annotations.secretome file and annotations.transmembrane file as I slightly altered the criteria for what's considered secretome / transmembrane which differs from Funannotate. I also added TMHMM, which Funannotate currently doesn't run.

2. Species_proteins.fasta 

3. Orthofinder output:

      * Orthogroups.txt
      
          Should look like this: 
      
          OG000000: Spec1_0001-T1 Spec2_4954-T1 Spec3_0492-T1   
          OG000001: Spec2_5928-T1 Spec3_1134-T1 Spec4_0031-T1   
          OG000002: Spec1_0353-T1 Spec3_3923-T1 Spec4_2953-T1   
          
          The script follows protein_id / locus_tag rules from NCBI (i.e. Uniq_####) where 'Uniq' is a unique identifier for the species and the #### following the underscore are the ID for the protein (if followed by -T or -mRNA) or gene (if using locus tags). Just be sure that the ID's used in Orthofinder (which will be the FASTA headers used in the Orthofinder run) match the first columns of the Funannotate annotations files so that proteins can be matched to the annotations files.

## Input file structure

```
Parent_directory
    |
    |
    +--Orthogroups.txt
    |
    |
    +--Isolate1 
    |       |
    |       |
    |       +-- Isolate1_proteins.fasta
    |       |
    |       +-- (All annotation files) <- make sure analysis names are in filename (such as 'pfam', 'dbcan', 'antismash', etc)
    |
    |
    +--Isolate2
            |
            |
            +-- Isolate2_proteins.fasta
            |
            +-- (All annotation files)
```


