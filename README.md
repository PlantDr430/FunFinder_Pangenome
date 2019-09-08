# FunFinder_Pangenome
#### Data analysis of Funannotate and Orthofinder outputs to Pangenome results

## Summary

This script was written to produce many of the results and figures in Wyka et al. (in preparation). It takes output files from [Funannotate](https://github.com/nextgenusfs/funannotate) and the Orthogroups.txt output from [Orthofinder](https://github.com/davidemms/OrthoFinder) and conducts pangenomic analyses and produces results in the form of text files and figures.

Results will consist of:

Pangenome stats
![Pangenome stats](https://github.com/PlantDr430/images/blob/master/Pangenome_stats.png)

Pangenome curve - with power law regression fit
![Pangenome curve - with power law regression fit](https://github.com/PlantDr430/images/blob/master/Pangenome_curve.png)

Pangenome fluidity - computed with bootstraps (resampling)
![Pangenome fluidity - computed with bootstraps (resampling)](https://github.com/PlantDr430/images/blob/master/Pangenome_fluidity.png)

Protein length boxplots for pangenome category (lettering for significance with Kruskal-Wallis Test)
![Protein lengths](https://github.com/PlantDr430/images/blob/master/Protein_lengths.png)

Bar charts for all analyses respective showing proportion of clusters in each pangenomic category if >= 50% of species represented in the cluster have at least 1 protein with the respective annotation (% identify can be altered) (letters for significance with Fischer Exact test)
![Bar charts](https://github.com/PlantDr430/images/blob/master/secretome_pangenome_bar.png)

## Output files

All figures as well as text files associated with images for reproduction purposes. Within in isolate directory in the species_results directory there will be fasta files containing proteins from each analysis provided, as well as, overall_results file for each analysis. 

## Dependencies 

1. [Python 3](https://www.python.org/downloads/)
2. [EffectorP_2.0](http://effectorp.csiro.au/software.html)

## Usage

```
usage: ./Orthofinder_2_pangenome.py [options] -d directory -a analyses -o output

optional arguments:
  -h, --help                      show this help message and exit
  -d , --directory                Directory containing folders for processing
  -a  [ ...], --analyses  [ ...]  Analyses to compare separated by commas [tmhmm,phobius,signalp,effectors,merops,dbcan,pfam,iprscan,antismash]
  -o , --out                      Output folder
  -f_p , --fischer_pvalue         Pvalue cutoff for Fischer exact test [default: 0.05]
  -k_p , --kruskal_pvalue         Pvalue cutoff for Kruskal-Wallis test [default: 0.05]
  -p , --percent                  Percent of isolates within a cluster that contains proteins with analysis hit,to consider the cluster as such [default: 0.50]
  -b , --bootstrap                Bootsraps for genome fluidity calculations [default: 5]
```

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
    +--Isolate  
    |       |
    |       |
    |       +-- Species_proteins.fasta
    |       |
    |       +-- (All annotation files) <- make sure analysis names are in filename (such as 'pfam', 'iprscan', 'dbcan', etc)
    |
    |
    +--Isolate
            |
            |
            +-- Species_proteins.fasta
            |
            +-- (All annotation files)
```

