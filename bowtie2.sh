#!/bin/bash
myname=`basename $0`
DEBUG=0
#DEFAULT_SAMPLEID=NA
DEFAULT_THREADS=1

usage() { echo "###########################################################################################################
                   $myname: map paired end reads on a genomic sequence with bowtie2

Usage: $myname -r <reads_R1.fastq(.gz)> -s <reads_R2.fastq(.gz)> -g <genome.fa> -o <outputdir> [-t <ncpus>]
       
Options:
  -t INT number of threads [1]

Output (in outputdir): reads.bam (sorted), reads.bam.flagstat, reads.bam.bwtlog, reads.q10.bam (mapped and MAPQ>=10), reads.q10.bam.flagstat
WARNING: needs the genomic sequence to be previously indexed 
WARNING: needs the first and second read files to contain \"_R1\" and \"_R2\" respectively
###########################################################################################################" 1>&2; exit 1;}

#a=$DEFAULT_SAMPLEID
t=$DEFAULT_THREADS

while getopts ":r:s:g:o:t:" opt; do
    case "${opt}" in
#        s)
#	s=${OPTARG}
#	;;
        r)
	r=${OPTARG}
	;;
        s)
	s=${OPTARG}
	;;
        g)
	g=${OPTARG}
	;;
        o)
	o=${OPTARG}
	;;
        t)
	t=${OPTARG}
	;;
#        *)
#	usage
#	;;
    esac
done

if [ $DEBUG == 1 ];
then
#   echo " s = ${s} " 
   echo " r = ${r} " 
   echo " s = ${s} " 
   echo " g = ${g} "
   echo " o = ${o} "
   echo " t = ${t} "
fi

#fastq1=$s
#fastq2="${fastq1/_R1/_R2}" # reads (second in pair)
fastq1=$r
fastq2=$s
name=`basename $fastq1 | sed s/_R1.*//`
genome=$g
outdir=$o
sample=$a
ncpu=$t
#readgroup="@RG\tID:"$sample"\tSM:"$sample
tmpid=`date +%m%d%H%M%N` # to avoid conflicting temporary files
tmpoutdir=$outdir/tmp$tmpid
unsorted=$tmpoutdir/$name.unsorted.bam
sam=$tmpoutdir/$name.sam
bam=$outdir/$name.bam
mono=$tmpoutdir/$name.mono.bam
q10=$outdir/$name.q10.bam

# sam=$outdir/`basename $fastq1 %_R1.*}.unsorted.sam`
# bam=${sam%sam}bam
# sorted=${bam%.unsorted.bam}
log=$bam.bwtlog

if [[ ! -r "$fastq1" ]]
    then
    echo "ERROR: could not read input reads1 \"$fastq1\""
    usage
fi

if [[ ! -r "$fastq2" ]]
    then
    echo "ERROR: could not read input reads2 \"$fastq2\""
    usage
fi

if [[ ! -r "$genome" ]]
    then
    echo "ERROR: could not read genomic reference \"$genome\""
    usage
fi

if [[ ! "$outdir" ]] 
    then
    echo "ERROR: please provide a name for the output directory"
    usage
fi

echo "##### Starting $0 - `date`"

echo "## Creating output directory..."
mkdir -p $outdir
if [[ ! -d "$outdir" ]]
    then
    echo "ERROR: could not create the output directory \"$outdir\""
    usage
fi
echo "## ...done in $outdir" 

# needs the genome to be mapped
if [[ ! -r "$genome".1.bt2 ]] || [[ ! -r "$genome".2.bt2 ]] || [[ ! -r "$genome".rev.1.bt2 ]] || [[ ! -r "$genome".rev.2.bt2 ]]
    then
    echo "ERROR: could not read one of the index files from $genome.1.bt2 to $genome.rev.2.bt2"
    echo "Please create an index for the genomic reference before mapping"
    exit 2
fi

## safety check in case the genomic sequence has been updated but not the index (critical risk)
if [[ `date -r $genome +%s` -le `date -r $genome.1.bt2 +%s` ]]
    then
    echo "## Using existing index file $genome.bwt..."
else
    echo "ERROR: index file \"$genome.1.bt2\""
    echo "       is older than the genome file \"$genome\""
    echo "Please regenerate the index or, if you are sure that the index corresponds to the genomic sequence, just run \"touch $genome.bwt\""
    exit 3
fi

echo "## Creating and moving to temporary directory..."
mkdir -p $tmpoutdir
if [[ ! -d "$tmpoutdir" ]]
    then
    echo "ERROR: could not create the temporary directory \"$tmpoutdir\""
    usage
fi
cd $tmpoutdir
echo "## ...done in $tmpoutdir" 

echo "## mapping with bowtie2 and converting to bam..."
#bwa sampe -s -n 20 $genome $sai1 $sai2 $fastq1 $fastq2 > $sam
#bwa mem -M -R "$readgroup" $genome $fastq1 $fastq2 1> $sam 2> $log
#bwa mem -M -R "$readgroup" -t $t $genome $fastq1 $fastq2 1> $sam 2> $log
bowtie2 -t -p $t -X 2000 -x $genome -1 $fastq1 -2 $fastq2 --met 60 --met-file $log -S $sam
samtools view $sam -Sbh -@ $t -o $unsorted
echo "## ...done in $unsorted"

echo "## sorting the reads..."
samtools sort -@ $t -m 4G $unsorted -O bam -o $bam
echo "## ...done in $bam"
#rm -f $unsorted

echo "## indexing the bam file..."
samtools index $bam
echo "## ...done"

echo "## computing summary stats on bam file..."
samtools flagstat $bam > $bam.flagstat
echo "## ...done in $bam.flagstat"

echo "## extracting mapped reads with MAPQ>=10 and mapped mate..."
samtools view -@ $t -F 12 -q 10 -bh $bam > $mono
echo "## ...done in $mono"

echo "## sorting reads by name to keep pairs with both reads mapped with mapq10 on the same chrom"
samtools sort -n -@ $t -m 4G $mono -o - | samtools view -h - | awk '/^@/{print;next} $7=="=" && $1==id {print prev;print $0}$7=="="{id=$1;prev=$0}' | samtools view -Sbh - > $unsorted.2
echo "## ...done in $unsorted"

echo "## sorting the reads..."
samtools sort -@ $t -m 4G $unsorted.2 -O bam -o $q10
echo "## ...done in $q10"
#rm -f $unsorted

echo "## indexing and computing summary stats on final bam file..."
samtools index $q10
samtools flagstat $q10 > $q10.flagstat
samtools idxstats $q10 > $q10.idxstats
echo "## ...done in $q10.flagstat"

# GET RID OF TEMPORARY FILES
rm -R $tmpoutdir

echo "##### Mapping done in $outdir - check $log and $q10.flagstat"
echo "##### End of $0 - `date`"
