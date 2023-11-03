# Metagenomics

Metagenomics workflow made as a part of an experiment designed to study the tempo and mode of adaptation of a new strain of *E. coli*, and other species present in the gut of obese mice. Those mice were leptin deficient mice and were colonized with a strain of *E. coli*. In addition to the invader strain, the mice were also colonized with another *E. coli* strain which typically resides in their normal gut microbiota.

The aim of the workflow is to process the metagenomics data and detect SNVs in metagenomic assembled genomes (MAGs).

The workflow was divided in three main bash scripts, that allow user's input to change certain parameters in the Assembly and inStrain pipeline.
<br/>
<br/>
<br/>
**Assembly.sh**

1)Pre-processing of fastq files

2)Remove host contamination

3)Run nonPareil

4)Co-assembly with Megahit

5)Quast

6)Read recruitment
<br/>
<br/>
<br/>
**binning.sh**

1)Create metadata for the SqueezeMeta pipeline

2)Binning using SqueezeMeta 
<br/>
<br/>
<br/>
**instrain.sh**

1)Run checkM

2)Changing checkM output file structure for dREP

3)Dereplication with dRep

4)Change contigs IDs with seqkit

5)ORFs prediction with prodigal

6)Align reads to MAGs with bowtie2

7)inStrain profile

8)Gene annotations with prokka

9)Check CDS in common from prokka and prodigal

10)Extract fom prokka's .gbk file the common CDS and information

11)Merge SNVs table with gene and gene product names

12)Merge taxonomic information from the MAGs with the SNVs table
<br/>
<br/>
<br/>
**Third party softwares used:**

megahit:
https://academic.oup.com/bioinformatics/article/31/10/1674/177884?login=true

NonPareil:
https://academic.oup.com/bioinformatics/article/30/5/629/246142?login=false

bbmap:
https://jgi.doe.gov/data-and-tools/software-tools/bbtools/

bowtie2:
https://www.nature.com/articles/nmeth.1923

SqueezeMeta:
https://www.frontiersin.org/articles/10.3389/fmicb.2018.03349/full

inStrain:
https://www.nature.com/articles/s41587-020-00797-0

dRep:
https://www.nature.com/articles/ismej2017126

sekit:
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163962

prokka:
https://academic.oup.com/bioinformatics/article/30/14/2068/2390517

prodigal:
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119




