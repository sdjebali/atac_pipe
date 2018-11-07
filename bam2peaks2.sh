#!/bin/bash
# myname=`basename $0` # this gave an error # ?? not in genotoul nor on the cluster
myname=$(basename "$0")
echo "##### `date` ####"
echo "# Starting command: $0 $@"

module load bioinfo/Java8

default_picard="/work/project/fragencode/tools/multi/picard-tools-2.1.1/picard.jar"
default_mt="MT"

#defining colours for logfile
red() { echo "$(tput setaf 1)$*$(tput setaf 9)"; }
green() { echo "$(tput setaf 2)$*$(tput setaf 9)"; }
blue() { echo "$(tput setaf 4)$*$(tput setaf 9)"; }

usage() { echo "###########################################################################################################

This script will: 
take a bam file, remove the mitochonrial reads, 
index the new bam file and indicate that the MT reads are gone by adding \"noMT\" before the extension,
calculate information using samtools idxstats and flagstat, & write this information to appropriately named files
remove duplicates using picard tools and create a new file with the nodup added before the extension,
index the new noMT.nodup.bam file,
calculate information using samtools idxstats and flagstat, & write this information to appropriately named files
call peaks using macs2, calculate multiple metrics using picard-tools.

Note that the grep arguement \"MT\" needs to be adapted to suit the nomenclature of the the genome being analysed.  e.g. for human hg19 it's \"chrM\"

Usage: $myname -b <mapped_reads.bam> -g <genome.fa> -o <outputdir> [-p <path/to/picard-tools> -m <default_mt_chrom_name>]
       
Options:
 -p : path to picard tools (default: $default_picard )
 -m : name of the mitochondrial chromosome (default: $default_mt )

Output (in outputdir): mapped_reads.noMT.bam, mapped_reads.noMT.nodup.bam, peaks/<macs_output>, summary metrics
CAUTION: please make sure the genome file was the one used for the mapping!
###########################################################################################################" 1>&2; exit 1;}

# usage example on the cluster:
#dir=/work2/project/fragencode ; qsub -N b2p_gt45 -cwd -q workq -j y -l mem=8G -o ~/work/tmp.qsub.b2p_gt45.`date +%Y%m%d%H%M%S`.log -b y "$dir/tools/bam2peaks.sh -b $dir/results/atacseq/capra_hircus/cd4/ATAC45_TAGGCATG_L004.q10.bam -g $dir/data/species/capra_hircus/CHIR_1.0.101/capra_hircus.fa -o $dir/results/atacseq/capra_hircus/cd4/ATAC45noMT -m 'gi|316926505|ref|NC_005044.2|'"

p=$default_picard
m=$default_mt

while getopts ":b:g:o:p:m:" opt; do
    case "${opt}" in
        b)
	b=${OPTARG}
	;;
        g)
	g=${OPTARG}
	;;
        o)
	o=${OPTARG}
	;;
        p)
	p=${OPTARG}
	;;
        m)
	m=${OPTARG}
	;;
    esac
done

input=$b
genome=$g
name=`basename $input | sed 's/.bam//'`
outdir=$o
picard=$p
mt=$m
peaksdir=peaks

nomt=$outdir/$name.noMT.bam
nodup=$outdir/$name.noMT.nodup.bam
peaks=$outdir/$peaksdir/$name.noMT.nodup.macs

log=$outdir/$myname.log

if [[ ! -r "$input" ]]
    then
    red "ERROR: could not read input bam file \"$input\""
    usage
fi

if [[ ! -r "$genome" ]]
    then
    red "ERROR: could not read genomic sequence \"$genome\""
    usage
fi

if [[ ! "$outdir" ]] 
    then
    red "ERROR: please provide a name for the output directory"
    usage
fi

blue "# Creating output directory..."
mkdir -p $outdir
if [[ ! -d "$outdir" ]]
    then
    red "ERROR: could not create the output directory \"$outdir\""
    usage
fi

cd $outdir
green "## ... moved to $outdir" 

blue "# ...Removing Mitochondrial reads"

# DEBUG
# echo "Running command: samtools idxstats $input | cut -f 1 | grep -v $mt | xargs samtools view -bh $input > $nomt"
samtools idxstats $input | awk '{print $1}' | grep -v $mt | xargs samtools view -bh $input > $nomt
green "## ...removing MT reads done in $nomt"

blue "#  ...Building nomt index"
#java -jar -Xmx4g $picard/BuildBamIndex.jar I=$nomt O=$nomt.bai VERBOSITY=ERROR
samtools index $nomt
green "## ...index for noMT done"

blue "# ...Getting summary stats on noMT.bam"
samtools idxstats $nomt > $nomt.idxstats
samtools flagstat $nomt > $nomt.flagstat
green "## ... summary stats for $nomt done"

blue "# ...Removing duplicates"
module load bioinfo/Java8
java -jar -Xmx4g $picard MarkDuplicates I=$nomt O=$nodup M=$nodup.metrics.txt REMOVE_DUPLICATES=true 
green "## ...removing duplicates done in $nodup (log in $nodup.metrics.txt )"

blue "# ...Building nomt.nodup index"
#java -jar -Xmx4g $picard/BuildBamIndex.jar I=$nodup O=$nodup.bai VERBOSITY=ERROR
samtools index $nodup
green "## ... index for $nodup done"

blue "# ...Getting summary stats on noMT.nodup.bam"
samtools idxstats $nodup > $nodup.idxstats
samtools flagstat $nodup > $nodup.flagstat
green "## ... summary stats for $nodup done"

blue "# ...Collecting multiple metrics"
java -jar -Xmx4g $picard CollectMultipleMetrics I=$nodup O=$nodup.CMM R=$genome 
green " ...multiple metrics done for $nodup"

blue "# ...Calling peaks using MACS2"
mkdir -p $outdir/$peaksdir
if [[ ! -d "$peaksdir" ]]
    then
    red "ERROR: could not create the peak directory \"$outdir/$peaksdir\""
    exit 2
fi

macs2 callpeak --nomodel -f BAMPE -t $nodup -n $peaks --keep-dup all --verbose 0
#macs2 callpeak --nomodel -f BAMPE -t $name.noMT.nodup.bam -n $peaks --keep-dup all --verbose 0
green "## ... peak calling done"

echo "##### End of process for $input - `date`"
