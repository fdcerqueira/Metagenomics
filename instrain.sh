#!/bin/bash
###directory of the bash script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

##metasqueeze project directory and sorted bams directory
project=$(cat /var/tmp/project.$PPID)

##read the output directory chosen during the assembly.sh script
START_DIR=$(cat /var/tmp/outputDir.$PPID)

##project name
project_name=$(cat /var/tmp/project_name.$PPID)
instrain=$START_DIR/instrain

#helper scripts
#format checkM output for dREP taxonomy_wf
format=$SCRIPT_DIR/format.py
#format checkM output for dREP lineage_wf
format2=$SCRIPT_DIR/format2.py

#create intermediate table from prokka .gbk file to add later the gene names and product
snvs=$SCRIPT_DIR/convert.py
##create permutations of strings for instrainc compare
#permutation=$SCRIPT_DIR/permutation.py
#change gbk file from prokka just to contain CDS in common with prodigal(instrain)
common_seq=$SCRIPT_DIR/common_seq.py
#take taxonomy from bins files and add to snvs instrain table just from assembly bins
taxonomy=$SCRIPT_DIR/update_scaffolds_bins.py
#take taxonomy from bins files and add to snvs instrain table just from assembly bins+database
taxonomy_all=$SCRIPT_DIR/update_scaffolds_all.py
#add gene names and products from prokka to snvs table
add_functions=$SCRIPT_DIR/add_function.R

function show_help() {

    echo "The script it will activate Instrain"
    echo "How to use: $0 parameters"
    echo ""
    echo "OPTIONS:"
    echo ""
    echo " -h < show help >"
    echo " -n < number of threads to be used >"
    echo " -p < number of processes to be used for instrain (default=6) >"
    echo " -ani < choose average nucleotide identity (ANI) for dereplication (e.g. 0.95 to 0.98) >"
    echo " -r_ani < read ANI, beteen 0-1 (default: 0.95). Read ANI should be at least 3% lower (-ani) than the expected differences between your reads and the reference genome >"
    echo " -min_cov < minimum coverage to call a mutation at inStrain profile (default: 5) >"
    echo " -min_scaf < Minimum number of reads mapping to a scaffold to proceed with profiling it (default: 5) >"
    echo " -min_gen >-Minimum number of reads mapping to a genome to proceed with profiling it (default: 0) >"
    echo ""
    echo "from the options below choose (if) you want $0 to perform them:"
    echo ""
    echo " -all < create instrain database with your bins and iMGMC - integrated Mouse Gut Metagenomic Catalog (https://github.com/tillrobin/iMGMC) >"
    echo " -mag < only use the co-assembly as database for inStrain >"
    echo " -lineage < use the lineage_wf modudule from checkM instead of taxonomy_wf, more accurate but much slower >"
    echo ""
    echo " e.g. -n 20 -ani 0.98 -r_ani 0.95 -min_cov 10 -min_scaf 10 -min_gen 5 -all"

return 0
}

###input commands
while [ ! -z "$1" ]
 do
      case "$1" in
      -n) shift 1;threads=$1 ;;
      -p) shift 1;proce=$1;;
      -ani) shift 1; ani=$1;;
      -r_ani)shift 1; r_ani=$1;;
      -min_cov)shift 1; min_cov=$1;;
      -min_scaf)shift 1; min_scaf=$1;;
      -min_gen)shift 1; min_gen=$1;;
      -h) shift 1; show_help; exit 1 ;;
      -all) shift 2; an=true;;
      *) break;;
    esac
shift
done
#lack of input errors
if [ -z "$threads" ] || [ -z "$ani" ] #|| [ -z "$r_ani" ] || [ -z "$min_cov" ] || [ -z "$min_scaf" ]
then
 show_help
 exit 1
fi

##overwrite output dir
if [ ! -d $instrain ]
 then
 	mkdir $instrain
 else
    	echo "output directory  exist. Do you want overwrite it? (y/n)"
    	read sim
        	if [ $sim = "n" ]
            	then
            	exit 1
    	   else
        	rm -rf $instrain
        	echo "Creating/overwriting ${instrain} directory"
        	mkdir $instrain
        fi
fi

###E.coli fasta
if [ -f "$SCRIPT_DIR/4_residentBacterialGenome.fa" ] && [ -f "$SCRIPT_DIR/1_invaderBacterialGenome.fa" ]
  then
    resident=$SCRIPT_DIR/4_residentBacterialGenome.fa
    invader=$SCRIPT_DIR/1_invaderBacterialGenome.fa
fi

#directories
output_bins=$instrain/databases
merged_genomes_set=$output_bins/MergedGenomeSet
dereplicated=$output_bins/dereplicated
bt2=$output_bins/bt2
sams=$output_bins/sams
profile=$instrain/profile
#compare=$instrain/compare
host_removal=${START_DIR}/host_removal
genomes_checkm=${output_bins}/genomes_checkm
tmp_dir_check=${output_bins}/tmp_checkm
checkm=$output_bins/checkM
checkm_final=$checkm/final
prokka=$merged_genomes_set/prokka
prodigal=$dereplicated/prodigal
intermediate_tables=$profile/intermediate_tables_SNPs
final_tables_SNPs=$profile/final_tables_SNPs
intermediate_function=$instrain/intermediate_function

#create directories
mkdir $output_bins
mkdir $dereplicated
mkdir $merged_genomes_set
mkdir $bt2
mkdir $sams
mkdir $genomes_checkm
mkdir $tmp_dir_check
mkdir $checkm
mkdir $checkm_final
mkdir $profile
#mkdir $compare
mkdir $prodigal
mkdir $final_tables_SNPs
mkdir $intermediate_tables
mkdir $intermediate_function

echo ""
echo "Do you have inStrain, bowtie2,samtools, prokka, seqkit, fastANI and openssl=1.0 installed in a conda env? (y/n)?"
read CONT
echo ""
if [ "$CONT" = "y" ]
then
  echo "Proceeding with the pipeline"
else
###choose to create conda env
### prompt for conda env name
echo ""
echo "Creating new conda environment, choose name"
read input
echo "Name ${input} was chosen"

#install fastani manually
wget -O $SCRIPT_DIR/fastANI-Linux64-v1.33.zip "https://github.com/ParBLiSS/FastANI/releases/download/v1.32/fastANI-Linux64-v1.33.zip"
unzip $SCRIPT_DIR/fastANI-Linux64-v1.33.zip
echo "export PATH=$SCRIPT_DIR/:\$PATH" >> ~/.bashrc #fastani path to .bashrc
source ~/.bashrc

#list name of packages
echo "installing base packages"
eval "$(conda shell.bash hook)"
conda create -y --name $input python=3.9
#search for conda to be able to activate it and allow to install the required programs
conda activate $input
#talvez pelo git
pip install instrain
conda install -y -c bioconda drep
conda install -y -c bioconda samtools=1.9 --force-reinstall
pip uninstall biopython
pip install biopython==1.74

#seqkit env
conda create -y -n seqkit
conda activate seqkit
conda install -y -c bioconda seqkit

#bowtie env
conda create -y -n bowtie2
conda activate bowtie2
conda install -y -c bioconda bowtie2

#manual change this conda shit
conda create -y -n checkm python=3.9
conda activate checkm
conda install -y numpy matplotlib pysam
conda install -y hmmer prodigal pplacer
pip install checkm-genome

if [ ! -d "$SCRIPT_DIR/checkm_data" ]
  then
    mkdir checkm_data
    wget -O $SCRIPT_DIR/checkm_data/checkm_data_2015_01_16.tar.gz "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
    checkm data setRoot $SCRIPT_DIR/checkm_data
fi

#env for prokka
conda create -y -n prokka-env
conda activate prokka-env
conda install -c conda-forge -c bioconda -c defaults prokka
###remove conflict of biopython with inStrain


##samtools conflict (shared library libncurses5 not present)
#sudo apt-get install libncurses5

###install prokka via brew
#/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
#brew install brewsci/bio/prokka
##remove minced dependency, it usually gives errors
#brew remove minced
#if there is an error with R in conda, try to change the version installed (more recent)
echo "conda env ${input} created, and the required packages as well"
fi

###activate the conda environment and select the environment you have fastp, bbmap, megahit
eval "$(conda shell.bash hook)"
echo "the list of your conda environments"
conda env list
echo "Write the one you wanted"
read conda
conda activate $conda
echo "Conda environment ${conda} activated"

#build database with mags and public database
if [ "$an" = "true" ] && [ ! -d $SCRIPT_DIR/dereplicated_genomes/ ]
  then
    echo "The public database available is the :"
    echo "Integrated Mouse Gut Metagenomic Catalog (https://doi.org/10.1016/j.celrep.2020.02.036)"
    echo "Downloading MAGs gut mouse database to $SCRIPT_DIR "
        wget -O $SCRIPT_DIR/iMGMC-hqMAGs-dereplicated_genomes.tar.gz "https://zenodo.org/record/3631711/files/iMGMC-hqMAGs-dereplicated_genomes.tar.gz"
        wget -O $SCRIPT_DIR/MAG-annotation_CheckM_dRep_GTDB-Tk.tar.gz "https://zenodo.org/record/3631711/files/MAG-annotation_CheckM_dRep_GTDB-Tk.tar.gz"
        ##move to metasqueeze final bins directory
        ##untar general database
        tar xvzf $SCRIPT_DIR/MAG-annotation_CheckM_dRep_GTDB-Tk.tar.gz -C $SCRIPT_DIR
        tar xvzf $SCRIPT_DIR/iMGMC-hqMAGs-dereplicated_genomes.tar.gz -C $SCRIPT_DIR

    if [ $? -ne 0 ]
      then
        echo "Unable to download the databases, do it manually from: https://github.com/tillrobin/iMGMC"
        exit 1
    fi
fi

###copy final bins with decided completion and redundancy
echo "creating lists of the bins generated by squeezemeta"
touch $project/results/DAS/bins.txt
awk 'BEGIN {FS= "\t"} $12 >= 50 && $13 <= 10 {print $1} ' \
$project/results/DAS/${project_name}_DASTool_summary.txt > $project/results/DAS/bins.txt


#add path to list
awk '$1="'$project/results/DAS/${project_name}_DASTool_bins/'"$1' $project/results/DAS/bins.txt > $output_bins/bins.txt

##copy bins
echo "copying all bins"
for i in `cat $output_bins/bins.txt`
  do
    cp ${i}.contigs.fa $genomes_checkm
done

#copy resident and invader e.coli
if [ ! -z "$resident" ] && [ ! -z "$invader" ]
  then
    echo "copying E.coli genomes Inv/Res"
    cp $resident $genomes_checkm
    cp $invader  $genomes_checkm
fi

if [ "$an" = "true" ]
  then
    echo "copying database genomes"
    cp $SCRIPT_DIR/dereplicated_genomes/*.fa $genomes_checkm
fi

eval "$(conda shell.bash hook)"
conda activate checkm
###run checkM separatly. There is a conflict with dRep
if [ "$an" = "true" ]
  then
checkm taxonomy_wf \
        domain "Bacteria" \
        -t ${threads} \
        -x fa \
        --tmpdir $tmp_dir_check \
        $genomes_checkm \
        ${checkm}
  else
    checkm lineage_wf \
    -t ${threads} \
    -x fa \
    --tmpdir $tmp_dir_check \
    $genomes_checkm \
    ${checkm}
fi

###get table
checkm qa $checkm/*.ms ${checkm} \
--file ${checkm_final}/bins.csv \
--threads ${threads} \
-o 1


##python file
echo "changing checkM output file structure for dREP"
if [ "$an" = "true" ]
  then
    python $format ${checkm_final}/bins.csv > ${checkm_final}/checkm_out.csv
  else
    python $format2 ${checkm_final}/bins.csv > ${checkm_final}/checkm_out.csv
fi

eval "$(conda shell.bash hook)"
conda activate instrain

#dereplicate bins for large amount of genomes
if [ "$an" = "true" ]
  then
    echo "dereplicating bins with greedy algorithms"
    dRep dereplicate $dereplicated \
    -g ${genomes_checkm}/*.fa \
    --S_algorithm fastANI \
    --multiround_primary_clustering \
    --greedy_secondary_clustering \
    --run_tertiary_clustering \
    --genomeInfo ${checkm_final}/checkm_out.csv \
    -comp 50 \
    -con 10 \
    -ms 10000 \
    -pa 0.9 \
    -sa ${ani} \
    -nc 0.30 \
    -cm larger \
    -p ${threads}
else
  echo "dereplicating co-assembly MAGs"
  dRep dereplicate $dereplicated \
  -g ${genomes_checkm}/*.fa \
  --S_algorithm ANImf \
  --genomeInfo ${checkm_final}/checkm_out.csv \
  -comp 50 \
  -con 10 \
  -ms 10000 \
  -pa 0.9 \
  -sa ${ani} \
  -nc 0.30 \
  -cm larger \
  -p ${threads}

fi

eval "$(conda shell.bash hook)"
conda activate seqkit
#rename headers of sequences .fa files (they are too long to be uploaded and be processed later with seqIO)
#rename headers of sequences .fa files (they are too long to be uploaded and be processed later with seqIO)

for i in $dereplicated/dereplicated_genomes/maxbin*.fa
    do
      filename=$(basename $i)
      echo "Changing sequence IDs of $filename"
      filename1=${filename%%.contigs.fa}
      filename2=`echo ${filename1} | grep -o -E "[0-9]+"`
      seqkit replace -p ".+" -r "${filename2}_{nr}" $i > $dereplicated/${filename1}.rename.fa
        #seqkit rmdup -n $dereplicated/${filename1}.renam.fa -D $dereplicated/duplicates.fa -o $dereplicated/${filename1}.rename.fa
done

if [ $? -ne 0 ]
  then
    echo "There were no maxbin Mags at $dereplicated/dereplicated_genomes"
fi


for i in $dereplicated/dereplicated_genomes/meta*.fa
  do
    filename=$(basename $i)
    echo "Changing sequence IDs of $filename"
    filename1=${filename%%.contigs.fa}
    filename2=`echo ${filename1} | grep -o -E "[0-9]+"`
    seqkit replace -p ".+" -r "${filename2}_{nr}" $i > $dereplicated/${filename1}.rename.fa
done

if [ $? -ne 0 ]
  then
  echo "There were no metabat mags $dereplicated/dereplicated_genomes"
fi


if [ "$an" = "true" ]
  then
    for i in $dereplicated/dereplicated_genomes/extra*.fa
      do
        filename=$(basename $i)
        echo "Changing sequence IDs of $filename"
        filename1=${filename%%.fa}
        filename2=`echo ${filename1} | tail -c5`
        seqkit replace -p ".+" -r "${filename2}_{nr}" $i > $dereplicated/${filename1}.rename.fa
    done

    for i in $dereplicated/dereplicated_genomes/iMGMC*.fa
      do
        filename=$(basename $i)
        echo "Changing sequence IDs of $filename"
        filename1=${filename%%.fa}
        filename2=`echo ${filename1} | tail -c5`
        seqkit replace -p ".+" -r "${filename2}_{nr}" $i > $dereplicated/${filename1}.rename.fa
    done

    for i in $dereplicated/dereplicated_genomes/single*.fa
      do
        filename=$(basename $i)
        echo "Changing sequence IDs of $filename"
        filename1=${filename%%.fa}
        filename2=`echo ${filename1} | tail -c5`
        seqkit replace -p ".+" -r "${filename2}_{nr}" $i > $dereplicated/${filename1}.rename.fa
    done
fi


filename=$(basename $dereplicated/dereplicated_genomes/4_residentBacterialGenome.fa)
filename1=${filename%%.fa}
filename2=`echo ${filename1} | tail -c5`
seqkit replace -p ".+" -r "${filename2}_{nr}" $dereplicated/dereplicated_genomes/4_residentBacterialGenome.fa > $dereplicated/${filename1}.rename.fa

if [ -s "$dereplicated/4_residentBacterialGenome.rename.fa" ]
  then
    rm $dereplicated/4_residentBacterialGenome.fa.rename.fa
fi

filename=$(basename $dereplicated/dereplicated_genomes/1_invaderBacterialGenome.fa)
filename1=${filename%%.fa}
filename2=`echo ${filename1} | tail -c5`
seqkit replace -p ".+" -r "${filename2}_{nr}" $dereplicated/dereplicated_genomes/1_invaderBacterialGenome.fa > $dereplicated/${filename1}.rename.fa

if [ -s "$dereplicated/1_invaderBacterialGenome.rename.fa" ]
  then
    rm $dereplicated/1_invaderBacterialGenome.rename.fa
fi

eval "$(conda shell.bash hook)"
conda activate instrain
##concatenate list of all bins
for i in $dereplicated/*.rename.fa
  do

    filename=$(basename $i)
    filename1=${filename%%.fa}

    prodigal -i $i -c -m -g 11 -p single -f sco -q -d $prodigal/${filename1}.fna
done


echo "merging .fna files"
cat $prodigal/*.fna > $merged_genomes_set/all_genomes.fna

echo "concatenating genomes"
cat $dereplicated/*.rename.fa > $merged_genomes_set/all_genomes.fa


##create horizontal list with all bins
echo "creating a list of the genomes"
ls ${dereplicated}/*.rename.fa | tr "\n" " " > ${dereplicated}/genomes.txt
lista=$(cat ${dereplicated}/genomes.txt)

#create stb file
echo "creating stb file"
parse_stb.py --reverse -f $lista -o $merged_genomes_set/genomes.stb

#create bowtie2 index
eval "$(conda shell.bash hook)"
conda activate bowtie2

echo "creating bowtie2 index"
if [ "$an" = "true" ]
  then
    bowtie2-build $merged_genomes_set/all_genomes.fa $output_bins/bt2/all_genomes --large-index --threads ${threads}
  else
    bowtie2-build $merged_genomes_set/all_genomes.fa $output_bins/bt2/all_genomes --threads ${threads}
fi


echo "mapping the reads to the genomes"

for i in $host_removal/cleaned*_R1.fastq.gz
  do
      filename=$(basename $i)
      filename1=${filename%%_R1.fastq.gz}
      filename2=${filename1%%_R1.fastq.gz}_R2.fastq.gz
      filename3=${filename1##cleaned-M_t-}

      bowtie2 -p ${threads} \
      -x $output_bins/bt2/all_genomes \
      -1 $host_removal/${filename} \
      -2 $host_removal/${filename2} \
      > $sams/${filename3}.sam
done

eval "$(conda shell.bash hook)"
conda activate instrain

###instrain profilling
echo "running Instrain profile"
for i in $sams/*.sam
  do
    filename=$(basename $i)
    filename1=${filename%%.sam}

###run database mode if bins+database
    if [ "$an" = "true" ]
      then
        inStrain profile \
        $i \
        $merged_genomes_set/all_genomes.fa \
        -o $profile/${filename1} \
        -p ${proce} \
        -g $merged_genomes_set/all_genomes.fna \
        -s $merged_genomes_set/genomes.stb \
        --database_mode \
        --min_cov ${min_cov} \
        --min_scaffold_reads ${min_scaf} \
        --min_genome_coverage ${min_gen} \
        --skip_plot_generation  2>&1 | tee -a $profile/$filename1.log
#run normal mode
      else
        inStrain profile \
        $i \
        $merged_genomes_set/all_genomes.fa \
        -o $profile/${filename1} \
        -p ${proce} \
        -g $merged_genomes_set/all_genomes.fna \
        -s $merged_genomes_set/genomes.stb \
        -l ${r_ani} \
        --min_cov ${min_cov} \
        --min_scaffold_reads ${min_scaf} \
        --min_genome_coverage ${min_gen} 2>&1 | tee -a $profile/$filename1.log
      fi

      if [ $? -ne 0 ]
        then
          echo "Error instrain"
          exit 1
      fi
done


##check for gzip Files of instrain profile output
for i in $sams/*.sam
  do
    filename=$(basename $i)
    filename1=${filename%%.sam}

    if [ -f "$profile/${filename1}/output/${filename1}_SNVs.tsv.gz" ]
      then
        gzip -dk $profile/${filename1}/output/${filename1}_SNVs.tsv.gz
    fi
done

###Add taxonomy of the bins and database bins in the final mutations table
echo "Adding taxonomy, gene names and products to the inStrain profile SNVs tables"
echo ""

echo "Step 1: gene names and products"

echo "Running prokka"
eval "$(conda shell.bash hook)"
conda activate prokka-env
#running proka for annotations
for i in $(find $dereplicated/ -type f -name "*.rename.fa")
  do

    filename=$(basename $i)
    filename1=${filename%%.fa}

    prokka --outdir $prokka --force --norrna --notrna --prefix $filename1 --cpus 0 $i
done

eval "$(conda shell.bash hook)"
conda activate instrain

##removing new annotations for resident and invader ecoli with the Instrain advised parameters for prodigal
#if [ -f "$prokka/4_residentBacterialGenome.rename.gbk" ] || [ -f "$prokka/1_invaderBacterialGenome.rename.gbk" ]
#  then
#    rm $prokka/4_residentBacterialGenome.rename.gbk
#    rm $prokka/1_invaderBacterialGenome.rename.gbk
#fi

echo "merging .gbk file of all bins"
cat $prokka/*.gbk > $merged_genomes_set/all_genomes.gbk
echo ""
#echo "Appending Resident and Invader E.coli gbk"
#if [ -s "$dereplicated/dereplicated_genomes/4_residentBacterialGenome.rename.fa" ]
#  then
#    cat $SCRIPT_DIR/ResidentGenomic.gbk >> $merged_genomes_set/all_genomes.gbk
#fi

#if [ -s "$dereplicated/dereplicated_genomes/4_residentBacterialGenome.rename.fa" ]
#  then
#    cat $SCRIPT_DIR/NC_000913.2.gbk >> $merged_genomes_set/all_genomes.gbk
#fi

echo "Merging .ffn files of all bins"
cat $prokka/*.ffn > $merged_genomes_set/all_genomes_p.ffn

echo "Checking if the CDS files from prodigal(inStrain) and prokka are equal"

prodi=$(grep ">" $merged_genomes_set/all_genomes.fna | wc -l)
prok=$(grep ">" $merged_genomes_set/all_genomes_p.ffn | wc -l)


eval "$(conda shell.bash hook)"
conda activate seqkit
echo "seqkit"
seqkit common -s $merged_genomes_set/all_genomes_p.ffn $merged_genomes_set/all_genomes.fna -o $merged_genomes_set/common.fasta 2>&1 | tee -a $merged_genomes_set/nr_seqs.txt
#&> $merged_genomes_set/nr_seqs.txt

##get the number of equal CDS sequences in both Files
sed -n 4p $merged_genomes_set/nr_seqs.txt | grep -o -E "[0-9]+" > $merged_genomes_set/nrr_seqs.txt

 p1=$(awk "FNR == 2 {print $1} " $merged_genomes_set/nrr_seqs.txt)
 p2=$(awk "FNR == 4 {print $1} " $merged_genomes_set/nrr_seqs.txt)

if [ "$p1" == "$p2" ]
  then
    seq_comp=$p1
  else
    seq_comp="not_equal"
fi

###python script to change the instrain mutations files: to add genes names, product, strand (+/-), and taxonomy
if [ "$p1" != "$prok" ] || [ "$p1" != "$seq_comp" ] || [ "$p2" != "$seq_comp" ]
  then
    echo "CDS sequences from prodigal and prokka are not the same."

    echo "the number of cds from prodigal: $prodi"
    echo "the number of CDS from prokka: $prok"
    echo "the number of equal CDSs are $p2 from a total of $p1"

fi

###change gbk file from prokka, in order to just contain the common CDS with prodigal(instrain)
echo " Step 1: Genes and gene products"
echo "extracting common CDS to .gbk"
eval "$(conda shell.bash hook)"
conda activate instrain

python $common_seq $merged_genomes_set/common.fasta $merged_genomes_set/all_genomes.gbk $merged_genomes_set/common.gbk
echo "creating intermediate table to latter add function and taxonomy"
python $snvs $merged_genomes_set/common.gbk $intermediate_function

#Rscript to merge snp table with gene, product names
for i in $profile/*/output/*_SNVs.tsv
  do

    filename=$(basename $i)
    filename1=${filename%%.SNVs.tsv}
    echo "Merging genes and products names to $filename1"
    Rscript --vanilla $add_functions $intermediate_function/snp.txt $i $final_tables_SNPs
done

echo " Step 2: Taxonomy"

for i in $final_tables_SNPs/*.csv
  do
    echo " Merging Taxonomy to the the file $i"

    if [ "$an" = "true" ]
      then
        python $taxonomy_all \
        $SCRIPT_DIR/gtdbtk.bac120.summary.tsv \
        $i \
        $project/results/DAS/${project_name}_DASTool_scaffolds2bin.txt \
        $project/results/tables/${project_name}.bin.tax.tsv \
        $merged_genomes_set/genomes.stb \
        $i \
        $profile
      else
        python $taxonomy \
         $project/results/tables/${project_name}.bin.tax.tsv \
         $project/results/DAS/${project_name}_DASTool_scaffolds2bin.txt \
         $i \
         $i \
         $final_tables_SNPs
    fi

done

echo ""
echo "instrain script ran successfully"
echo "exiting"
exit 0