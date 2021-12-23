#!/bin/bash
###directory of the bash script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

###function with the parameters
function show_help() {

    echo "The script it will activate conda and run quality control (fastp)\n"
    echo "remove host contamination (bbduk.sh), diversity coverage estimations(NonPareil)"
    echo "assembly (MEGAHIT),assembly quality control (metaquast)"
    echo "and read recruitment to the assembly (bbwrap.sh)"
    echo "How to use: $0 parameters"
    echo ""
    echo "OPTIONS:"
    echo ""
    echo " -h <show help>"
    echo " -n <number of threads to be used>"
    echo " -xmx < memory RAM to be used with bbduk and bbwrap | e.g -xmx 3g for 3 gigabytes"
    echo " -r < fraction of the computer's total memory [0-1] to be used by MEGAHIT"
    echo " -i <input folder with the raw fastq.gz files>"
    echo " -o <output directory>"
    echo " -k_min <kmer minimum size for MEGAHIT assembly, it must be odd number, the default is 21>"
    echo " -k_max <kmer maximum size for MEGAHIT assembly,  must be odd number the default is 141>"
    echo " -k_step <increment of kmer size of each iteration, must be even number "
    echo ""
    echo " -no-mercy < write this option in case to not add mercy kmers. IF YOU CHOOSE THIS, you still have to declare -k_min, -k_max, and -k_step"
    echo " -default  < option to run MEGAHIT with default parameters. IF YOU CHOOSE THIS THERE IS NO NEED TO WRITE THE OTHER MEGAHIT PARAMETERS"
    echo ""
    echo " e.g. ./assembly.sh  -i /path/to/fastq -o path/output -n 10 -r 0.8  -xmx 3g -default_megahit"
    echo " e.g. ./assembly.sh  -i /path/to/fastq -o path/output -n 10 -r 0.8  -xmx 3g -k_min 37 -k_max 129"
    echo " e.g. ./assembly.sh  -i /path/to/fastq -o path/output -n 10 -r 0.8  -xmx 3g -k_min 37 -k_max 129 -nomercy"
return 0
}

###input parameters
while [ ! -z "$1" ]
 do
      case "$1" in
      -i) shift 1;inputFolder=$1 ;;
      -o) shift 1;outputDir=$1 ;;
      -r)shift 1;  RAM=$1 ;;
      -n) shift 1; CORES=$1 ;;
      -xmx)shift 1;  mem=$1;;
      -k_min) shift 1; declare -i KMIN=$1;;
      -k_max) shift 1; declare -i KMAX=$1;;
      -k_step) shift 1; declare -i KSTEP=$1;;
      -h) shift 1; show_help; exit 1 ;;
      -nomercy)  shift 2; nomercy=true;;
      -default_megahit) shift 2; default=true;;
      *) break;;
    esac
shift
done

##check if base input is present
if [ -z "$inputFolder" ] || [ -z "$outputDir" ] || [ -z "$RAM" ] || [ -z "$CORES" ] || [ -z "$mem" ]
 then
	show_help
  exit 1
fi

##outputdir
if [ ! -d $outputDir ]
 then
 	mkdir $outputDir
 else
    	echo "output directory  exist. Do you want overwrite it? (y/n)"
    	read sim
        	if [ $sim = "n" ]
            	then
            	exit 0
    	else
        	rm -rf $outputDir
        	echo "Creating/overwriting ${outputDir} directory"
        	mkdir $outputDir
        fi
fi


#host reference and nonpareil Rscript
mouseREF=${SCRIPT_DIR}/GCA_003774525.2_ASM377452v2_genomic.fna
pareilR=${SCRIPT_DIR}/nonPareil.R

#declare paths for input and outputs
#raw fastq files and its folder
raw_fastq=${inputFolder}
##fastp input and output
fastp_output=${outputDir}/filtered_reads
fastp_reports=${outputDir}/fastp_reports
###merge
merged_reads=${fastp_output}/merged_reads
##bbsplit/bbduk output
host_removal=${outputDir}/host_removal
all_reads=${host_removal}/all_reads
bbduk_reports=${host_removal}/reports_bbduk
contaminated_reads=${host_removal}/condaminated_reads
#nonPareil
non_Pareil=${outputDir}/non_Pareil/
non_Pareil_reports=${non_Pareil}/NP_reports/
##MEgahit assembly
Assembly=${outputDir}/Assembly
co_assembly=${Assembly}/co-assembly
#quast reports
quast_reports=${outputDir}/quast_reports
###read recruitment
Assembly_alignments=${outputDir}/Assembly_alignments
reports_read_recruitment=${outputDir}/reports_read_recruitment
refs=${Assembly_alignments}/refs
sorted_bams=${Assembly_alignments}/sorted_bams
Assembly_read_recruitment=${outputDir}/Assembly_read_recruitment
mapped_reads=${Assembly_read_recruitment}/mapped_reads
unmapped_reads=${Assembly_read_recruitment}/unmapped_reads
binning=${outputDir}/Binning

##make all directories
mkdir  $fastp_reports
mkdir  $fastp_output
mkdir  $merged_reads
mkdir  $host_removal
mkdir  $bbduk_reports
mkdir  $contaminated_reads
mkdir  $non_Pareil
mkdir  $non_Pareil_reports
mkdir  $Assembly
mkdir  $co_assembly
mkdir  $quast_reports
mkdir  $Assembly_alignments
mkdir  $sorted_bams
mkdir  $refs
mkdir  $Assembly_read_recruitment
mkdir  $unmapped_reads
mkdir  $mapped_reads
mkdir  $reports_read_recruitment
mkdir  $all_reads
mkdir $binning
###exported variable to temporary files to be read by the other scripts
echo "$outputDir" > /var/tmp/outputDir.$PPID
echo "$sorted_bams" > /var/tmp/sorted_bams.$PPID

echo ""
echo "Do you have the required programs (fastp, bbmap, nonPareil, MEGAHIT, and samtools)  installed in a conda env? (y/n)?"
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
  echo "Name ${input} was chosen";
    #list name of packages
        echo "installing base packages"
        conda create -y --name $input python=3.6
#search for conda to be able to activate it and allow to install the required programs
        eval "$(conda shell.bash hook)"
        conda activate $input
        conda install -y -c bioconda/label/cf201901 bbmap
        conda install -y -c bioconda/label/cf201901 fastp
        conda install -y -c conda-forge r-base=4.1.0
        conda install -y -c bioconda/label/cf201901 megahit
        conda install -y -c bioconda samtools=0.1.19
        #conda env for nonpareil, to avoid conflicts
        conda create -y --name nonpareil python=3.6
        conda activate nonpareil
        conda install -y -c bioconda/label/cf201901 nonpareil
        #conda env for quast to avoid conflicts
        conda create -y --name quast python=3.6
        conda activate quast
        conda install -y -c bioconda/label/cf201901 quast


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

###the pipeline
echo ""
echo "Running fastp"
echo ""
##quality control
for f1 in $raw_fastq/*_R1_001.fastq.gz
do
    filename=$(basename $f1)
    filename1=${filename%%_R*_*}
    f2=${f1%%_R1_001.fastq.gz}_R2_001.fastq.gz
    fastp -q 20 \
    --length_required 100 \
    --detect_adapter_for_pe \
    -p 3 -5 \
    -M 20 \
    -W 4\
    -c \
    -i $f1 \
    -I $f2 \
    -o $fastp_output/t-$(basename $f1) \
    -O $fastp_output/t-$(basename $f2) \
    -j $fastp_reports/${filename1}.json \
    -h $fastp_reports/${filename1}.html
done

####depending on the file extension, it will merge samples or not

    echo ""
    echo "Merging reads from different lanes"
    echo ""
    for i in $(find $fastp_output -type f -name "*.fastq.gz" \
      | while read F; do basename $F \
      | rev \
      | cut -c 22- \
      | rev; done \
      | sort \
      | uniq)
        do
          echo "Merging ${filename} R1"
          filename=${i%%_L00*}
          cat $fastp_output/${filename}_L00*_R1_001.fastq.gz > $merged_reads/M_${filename}_R1.fastq.gz
          echo "Merging ${filename} R2"
          cat $fastp_output/${filename}_L00*_R2_001.fastq.gz > $merged_reads/M_${filename}_R2.fastq.gz
        done

###host_removal
echo ""
echo "Removing host contamination"
echo "results stats are printed to stderr"
echo "you can check them at: ${bbduk_reports}"

for f1 in $merged_reads/*_R1.fastq.gz
do
  filename=$(basename $f1)
  filename1=${filename%%_R*}
  f2=${f1%%_R1.fastq.gz}_R2.fastq.gz

  bbwrap.sh \
  -Xmx${mem} \
  t=$CORES \
  minid=0.9 \
  ref=${mouseREF} \
  in=$f1 \
  in2=$f2 \
  out=$host_removal/cleaned-$(basename $f1) \
  out2=$host_removal/cleaned-$(basename $f2) \
  outm=$contaminated_reads/mouse-$(basename $f1) \
  outm2=$contaminated_reads/mouse-$(basename $f2) \
  2> $bbduk_reports/${filename1}.txt
done

###nonPareil
echo ""
echo "Running nonPareil"
echo ""
echo "decompressing fastq files"

for i in $host_removal/*_R1.fastq.gz
  do

  filename=$(basename $i)
  filename1=${filename%%_R1.fastq.gz}

  zcat $i > $host_removal/${filename1}_R1.fastq
done

eval "$(conda shell.bash hook)"
conda activate nonpareil
##non pareil only uses R1 reads, not R2.
for i in $host_removal/*R1.fastq
do
  filename=$(basename $i)
  filename1=${filename%%_R*}
  nonpareil \
  -s $i \
  -T kmer \
  -f fastq \
  -t $CORES \
  -b $non_Pareil/$filename1
done

###run non pareil rscript for plots and reports
echo ""
echo "Creating nonPareil plots and tables"
echo "you can check them at: ${non_Pareil_reports}"

Rscript --vanilla ${pareilR} ${non_Pareil} ${non_Pareil_reports}

##check sucess of Rscript
if [ $? -eq 0 ]
then
  echo "nonPareil script executed"
else
  echo "Script exited with error." >&2
  echo ""
  echo "do you want to proceed without the nonPareil reports and diversity coverage estimations? (y/n)"
  read answer
  if [ $answer = "y" ]
    then
    echo "proceeding with assembly"
  else
    exit 1
  fi
fi


#activate assembly
eval "$(conda shell.bash hook)"
conda deactivate
conda activate $conda

###assembly
echo ""
echo "Assembly with MEGAHIT"
echo ""
echo "Merging Files for co-assembly"
###merge R1.fastq.gz
ls -1 $host_removal/*_R1.fastq.gz | \
while read fn
 do
  cat "$fn" >> $all_reads/all_reads_R1.fastq.gz
done
###merge R2.fastq.gz
ls -1 $host_removal/*_R2.fastq.gz | \
while read fn
 do
  cat "$fn" >> $all_reads/all_reads_R2.fastq.gz
done

###assembly
if [ "$default" = "true" ]
  then
    megahit \
    -1 $all_reads/all_reads_R1.fastq.gz \
    -2 $all_reads/all_reads_R2.fastq.gz \
    -o $Assembly/all_reads \
    -m $RAM \
    -t $CORES \
    --out-prefix co-assembly

elif [ "$nomercy" = "true" ]
  then
    megahit \
    -1 $all_reads/all_reads_R1.fastq.gz \
    -2 $all_reads/all_reads_R2.fastq.gz \
    -o $Assembly/all_reads \
    --k-min $KMIN \
    --k-max $KMAX \
    --k-step $KSTEP \
    -m $RAM \
    -t $CORES \
    --no-mercy \
    --min-contig-len 500 \
    --out-prefix co-assembly

else
    megahit \
    -1 $all_reads/all_reads_R1.fastq.gz \
    -2 $all_reads/all_reads_R2.fastq.gz \
    -o $Assembly/all_reads \
    --k-min $KMIN \
    --k-max $KMAX \
    --k-step $KSTEP \
    -m $RAM \
    -t $CORES \
    --min-contig-len 500 \
    --out-prefix co-assembly

    if [ ! -f "$Assembly/all_reads/*.fa" ]
      then
        echo "Assembly went wront. Exit"
        exit 1
    fi
fi

##activate env for quast
eval "$(conda shell.bash hook)"
conda deactivate
conda activate quast
#quast reports
echo ""
echo "Performing quast reports"
echo ""

for i in $Assembly/all_reads/*.contigs.fa
  do
    filename=$(basename $i)
    filename1=${filename%%.contigs.fa}
    metaquast ${i} -o $quast_reports/${filename1} -t $CORES
done

###bbwrap.sh to check read recruitment, and get sorted .bams. bbmap programs print to stderror
echo ""
echo "Evaluating read recruitment to the assembly, and generating sorted .bam"
echo "results stats are printed to stderr"
echo "you can check them at: ${Assembly_alignments}"

eval "$(conda shell.bash hook)"
conda deactivate
conda activate $conda

for i in $host_removal/cleaned*_R1.fastq.gz
  do
    filename=$(basename $i)
    filename1=${filename%%_R*}
    f2=${filename%%_R1.fastq.gz}_R2.fastq.gz
    filename2=${filename1##cleaned-M_t-}

    bbwrap.sh \
    -Xmx${mem} \
    t=$CORES \
    in=$i \
    in2=$host_removal/${f2} \
    ref=$Assembly/all_reads/co-assembly.contigs.fa \
    path=$refs/${filename1} \
    out=$sorted_bams/${filename2}.sam \
    outm=$mapped_reads/mapped-${filename} \
    outm2=$mapped_reads/mapped-${f2} \
    outu=$unmapped_reads/un-${filename1} \
    outu2=$unmapped_reads/un-${f2} \
    ambig=toss \
    2> $Assembly_alignments/${filename2}.txt \
    bamscript=bs.sh
    sh bs.sh
done

##make reports with % of reads mapped to assembly, using awk
echo ""
echo "Making report with % of reads mapped to assembly"
echo "you can check them at: ${reports_read_recruitment}"

for i in $Assembly_alignments/*.txt
  do
    filename=$(basename $i)
    filename1=${filename%%.txt}
    awk "/kBases/ {for(i=1; i<=44; i++) {getline; print}}" \
    $i > $reports_read_recruitment/${filename1}.txt
    awk "/mapped/ {getline; print} " \
    $reports_read_recruitment/${filename1}.txt > $reports_read_recruitment/i-${filename1}.txt
    awk  '{ sum += $2 } END { print "Reads mapped to assembly:" (sum / NR)"%" }' \
    $reports_read_recruitment/i-${filename1}.txt > $reports_read_recruitment/${filename1}.txt
done

echo ""
echo "cleaning intermediate files from folders"
#cleaned unmmerged reads
rm $fastp_output/*.fastq.gz
##host reads
rm -rf $contaminated_reads
##nonPareil
rm $host_removal/*.fastq
rm $non_Pareil/*.npl
rm $non_Pareil/*.npa
##read reports_read_recruitment
rm $reports_read_recruitment/i-*.txt

echo ""
echo "Assembly and reports done"
exit
