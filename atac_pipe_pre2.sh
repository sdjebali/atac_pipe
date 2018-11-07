# atac_pipe_pre2.sh
# difference with atac_pipe_pre.sh is that it uses bowtie2.sh and bam2peaks2.sh
# which themselves were changed because of the new version of bowtie2 (2.3.3.1 instead of 2.1.0)
# and because of the new path to picard tools and also maybe new way of specifying output
# for samtools sort and the fact that most bam files did not have the .bam extension with bowtie.sh

# Will make peaks from atacseq pe reads using bowtie and macs
# uses 4 threads for mapping and needs at least 40G of ram
# recommended to send it to cluster with 4 threads and 16G per thread
# a test sample is in 
# fastq1trim=/work/project/fragencode/workspace/sdjebali/fragencode/atacseq/pipeline/smalltest/oct2017/trimming/subset_1_ATAC52_atacseq_combined_R1_val_1.fq.gz
# fastq2trim=/work/project/fragencode/workspace/sdjebali/fragencode/atacseq/pipeline/smalltest/oct2017/trimming/subset_1_ATAC52_atacseq_combined_R2_val_2.fq.gz

# usage
# atac_pipe_pre2.sh fastq1 fastq genome outdir



# Check the input files and if OK assign to variables
#####################################################
if [ ! -e "$1" ] || [ ! -e "$2" ] || [ ! -e "$3" ] || [ ! -e "$4" ]
then
    echo "" >&2
    echo "          Make peaks from ATAC-seq PE reads using bowtie2 and macs2" >&2
    echo "" >&2
    echo USAGE: >&2
    echo "atac_pipe_pre2.sh fastq1 fastq genome outdir" >&2
    echo "" >&2
    echo "Takes as input:" >&2
    echo "- a fastq.gz file with atac-seq read1 named *_R1.fastq.gz" >&2
    echo "- a fastq.gz file with atac-seq read2 named *_R2.fastq.gz" >&2
    echo "- a fasta file of the genome that has been indexed with bowtie2" >&2
    echo "- an output directory" >&2
    echo "" >&2
    echo "Provides:" >&2
    echo "- several intermediate files as well as the peaks obtained from macs on the reads mapped with bowtie2 that were trimmed and to" >&2
    echo "  which the mitochondrial and duplicate reads were removed (q10 mappings)" >&2
    echo "" >&2
    echo "Important note:" >&2  
    echo "- cannot be run twice in the same directory without loosing previous outputs since uses fixed names for outputs" >&2
    echo "" >&2
    exit 1
else
    fastq1=$1
    fastq2=$2
    genome=$3
    outdir=$4
fi

# Variables from input files
############################
fastq1base=`basename $fastq1`
fastq2base=`basename $fastq2`

# Paths
#######
## general
###########
path="`dirname \"$0\"`"             # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
## binaries and scripts
#######################
binDir=$rootDir/bin
awkDir=$rootDir/src/awk
bashDir=$rootDir/src/bash
## Output file directories
###########################
# a. trimming
trimdir=$outdir/trimming
# b. fastqc
qcdir=$outdir/fastqc
# c. mapping
mapdir=$outdir/mapping
# d. peaks
peakdir=$outdir/peaks
# Programs and scripts
######################
# Bin
#####
trim=$binDir/trim_galore
fastqc=$binDir/fastqc
# Bash
######
map=$bashDir/bowtie2.sh
peak=$bashDir/bam2peaks2.sh



###########
## START  #
###########
header="Executing atac_pipe_pre.sh"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
pipelineStart=$(date +%s)

# Make a directory for the run if it does not exist as well as the necessary subdirs
####################################################################################
mkdir -p $outdir
# a. trimming
mkdir -p $trimdir
# b. fastqc
mkdir -p $qcdir
# c. mapping
mkdir -p $mapdir
# d. peaks
mkdir -p $peakdir


#################################
# Trim reads to remove adaptors #
#################################
# uses as input files called
# - ~/fragencode/data/reads/atacseq/sus_scrofa/liver/ATAC62/ATAC62_atacseq_combined_R1.fastq.gz
# - ~/fragencode/data/reads/atacseq/sus_scrofa/liver/ATAC62/ATAC62_atacseq_combined_R2.fastq.gz
echo I am trimming the reads to remove the adaptors... >&2
cd $trimdir
$trim --stringency 3 -q 20 --paired --nextera -o $trimdir $fastq1 $fastq2 > trim_galore.out 2> trim_galore.err
echo done >&2
# 2h30 with 1 cpu and 0.5 G used
# produces files called:
########################
# - ATAC62_atacseq_combined_R1_val_1.fq.gz
# - ATAC62_atacseq_combined_R2_val_2.fq.gz

fastq1trim=$trimdir/${fastq1base%.fastq.gz}_val_1.fq.gz
fastq2trim=$trimdir/${fastq2base%.fastq.gz}_val_2.fq.gz

######################################
# Check reads post-trim using fastqc #
######################################
echo I am checking the read quality with fastqc... >&2
cd $qcdir
$fastqc -o $qcdir --quiet --dir $qcdir $fastq1trim $fastq2trim > fastqc.out 2> fastqc.err
echo done >&2
# 16 minutes with 1 cpu and 1.6G of ram
# produces files called
# - ATAC62_atacseq_combined_R1_val_1.fq_fastqc.zip
# - ATAC62_atacseq_combined_R1_val_1.fq_fastqc.html
# - ATAC62_atacseq_combined_R2_val_2.fq_fastqc.zip
# - ATAC62_atacseq_combined_R2_val_2.fq_fastqc.html


#####################
# Map trimmed reads #
#####################
echo I am mapping the reads... >&2
cd $mapdir
$map -r $fastq1trim -s $fastq2trim -g $genome -o $mapdir -t 4 > bowtie.out 2> bowtie.err
echo done >&2
# 7 hours with 4 cpus and 25.5G of ram
# produces files called like this
# - ATAC62_atacseq_combined.bam.bwtlog
# - ATAC62_atacseq_combined.bam
# - ATAC62_atacseq_combined.bam.bai
# - ATAC62_atacseq_combined.bam.flagstat
# - ATAC62_atacseq_combined.q10.bam
# - ATAC62_atacseq_combined.q10.bam.bai
# - ATAC62_atacseq_combined.q10.bam.flagstat
# - ATAC62_atacseq_combined.q10.bam.idxstats


####################################################################
# Remove MT reads, remove duplicates, call peaks for q10.bam files #
####################################################################
# need to check whether .fastq.gz is needed as opposed to .fq.gz for trimming
# and also if we need _R1_val_1.fq.gz as input to bowtie
# another possibility is to specify the name of the bam in bowtie script?
echo I am removing the MT and duplicate reads and calling peaks from the q10 mappings... >&2
cd $peakdir
bam=$mapdir/${fastq1base%_R1.fastq.gz}.q10.bam
$peak -b $bam -g $genome -o $peakdir > bam2peaks.out 2> bam2peaks.err
echo done >&2
# 1 hour with 1 cpu and 7.5G
# produces files like this
# - ATAC62_atacseq_combined.q10.noMT.bam
# - ATAC62_atacseq_combined.q10.noMT.bam.bai
# - ATAC62_atacseq_combined.q10.noMT.bam.idxstats
# - ATAC62_atacseq_combined.q10.noMT.bam.flagstat
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.metrics.txt
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.bai
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.idxstats
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.flagstat
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.CMM.insert_size_metrics
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.CMM.insert_size_histogram.pdf
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.CMM.alignment_summary_metrics
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.CMM.quality_distribution_metrics
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.CMM.quality_distribution.pdf
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.CMM.quality_by_cycle_metrics
# - ATAC62_atacseq_combined.q10.noMT.nodup.bam.CMM.quality_by_cycle.pdf
# - peaks/ATAC62_atacseq_combined.q10.noMT.nodup.macs_peaks.xls
# - peaks/ATAC62_atacseq_combined.q10.noMT.nodup.macs_peaks.narrowPeak
# - peaks/ATAC62_atacseq_combined.q10.noMT.nodup.macs_summits.bed


##########
# Clean  #
##########
echo "I am cleaning" >&2
rm $mapdir/*.q10.bam
rm $mapdir/*.q10.bam.bai
rm $peakdir/*.q10.noMT.bam
rm $peakdir/*.q10.noMT.bam.bai
echo done >&2

########
# END  #
########
pipelineEnd=$(date +%s)
echo "atac_pipe_pre for $lid completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min " >&2

# disable extglob
shopt -u extglob

exit 0
