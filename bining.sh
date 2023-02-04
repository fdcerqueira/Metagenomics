#!/bin/bash

###directory of the bash script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

##read the output directory chosen during the assembly.sh scriptp
START_DIR=$(cat /var/tmp/outputDir.$PPID)

### Rfile to prepare the samples files for the metasueeze intput
sample_tests=${SCRIPT_DIR}/sample_tests.R
###direcotry
Assembly=${START_DIR}/Assembly
host_removal=${START_DIR}/host_removal
binning=${START_DIR}/Binning
##binning direcotry where everything will be outputed



function show_help() {

    echo "The script it will activate Metasqueeze and perform binning and functional analysis"
    echo "remove host contamination (bbduk.sh), diversity coverage estimations(NonPareil)"
    echo "How to use: $0 parameters"
    echo ""
    echo "OPTIONS:"
    echo ""
    echo " -h <show help>"
    echo " -n <number of threads to be used>"
    echo " -p < project name were the Metasqueeza files will be outputed, IT WILL CREATE A FOLDER WITH THE BINNING FOLDER ON YOUR ASSEMBLY PROJECT >"
    echo " -m <memory ram to be used by squeezemeta in Gb"
    echo " -o <output folder of the databases for squeezemeta. ONLY USE THIS OPTION IN CASE YOU NEED TO DOWNLOAD THE DATABASES (~ 200/250 GB)"
    echo ""
    echo " e.g. -n 4 -p meta -m 60 -o /path/to/database"
return 0
}

function samples() {
  echo ""
  echo "Confirm all the samples and names are given to the file test.samples."
  echo ""
  echo "If the file is not properly edited, go to the directory ${host_removal} and make the necessary changes"
  echo ""
  echo "it will wait until you write <y>"
}

while [ ! -z "$1" ]
 do
      case "$1" in
      -n) shift 1;threads=$1 ;;
      -p) shift 1;project=$1 ;;
      -m) shift 1; MEM=$1;;
      -h) shift 1; show_help; exit 1 ;;
      -o) shift 1; output_database=$1;;
      *) break;;
    esac
shift
done

if [ -z "$threads" ] || [ -z "$project" ] || [ -z "$MEM" ]
then
 show_help
 exit 1
fi

###variables to be stored at vat/tmp for the instrain script
echo "${binning}/${project}" > /var/tmp/project.$PPID
echo "${project}" > /var/tmp/project_name.$PPID

echo ""
echo "Do you have SqueezeMeta and compareM installed in a conda env? (y/n)?"
read CONT
echo ""
if [ "$CONT" = "y" ]
then
  echo "Proceeding with the pipeline"
else
###choose to create conda env
### prompt for conda env name
echo ""
echo "Creating new conda environment named SqueezeMeta"

        eval "$(conda shell.bash hook)"
        mamba create -y -n SqueezeMeta4 -c conda-forge -c bioconda -c fpusan squeezemeta=1.6
        mamba install -y -n SqueezeMeta4 -c bioconda comparem
	echo "conda env SqueezeMeta created, and the required packages as well"
fi
###activate the conda environment and select the environment you have
eval "$(conda shell.bash hook)"
conda activate SqueezeMeta7
echo "Conda environment ${conda} activated"

#activate perl from squeezemeta conda environment and not from the system  to avoid conflict between perl and cpan versions
echo "activating SqueezeMeta perl"
export PATH="${CONDA_PREFIX}/bin:$PATH"

##restar process in case of error from the last step of metasqueze, by detecting if the project folder exists
if [ ! -d $output_database ]
 then
  mkdir -p $output_database
  download_databases.pl ${output_database}
  echo "testing SqueezeMeta installation"
  test_install.pl
else
  if [ -z "$(ls -A ${output_database})" ]
    then
      download_databases.pl ${output_database}
      echo "testing SqueezeMeta installation"
      test_install.pl
    else
    	echo "The databases directory already exists and the databases are downloaded. Do you want overwrite it? (y/n)"
    	read sim
        	if [ $sim = "n" ]
            	then
            	echo "Proceding with the pipeline"
    	    else
        	   rm -rf $output_database
        	    echo "Creating/overwriting ${output_database} directory"
        	    mkdir -p $output_database
              download_databases.pl ${output_database}
              echo "testing SqueezeMeta installation"
              test_install.pl
          fi
    fi
fi


##listing all cleaned fastq to later create test.samples file
echo "creating sample.test file"
touch ${host_removal}/the.samples
ls ${host_removal}/*.fastq.gz > ${host_removal}/the.samples

##while loop to remove directory path from samples basename
while IFS= read -r line
do
    basename "$line"
done < ${host_removal}/the.samples > ${host_removal}/testt.samples

Rscript --vanilla ${sample_tests} ${host_removal}

#show samples on screen
cat ${host_removal}/test.samples
#show function samples
samples

read yes
if [ "$yes" = "y" ]
  then echo ""
else
  samples
  exit 1
fi

###metasqueeze workflow
echo "Proceeding with squeezemeta"
SqueezeMeta.pl -m coassembly \
-s $host_removal/test.samples \
-t ${threads} \
-canumem ${MEM} \
-p $binning/${project} \
-f ${host_removal} \
-extassembly ${Assembly}/all_reads/co-assembly.contigs.fa 


##generate tabular results to import to R

sqm2tables.py ${binning}/${project} ${binning}/${project}/results/tables
