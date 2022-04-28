

file=$1


f=$(echo $file | rev | cut -f1 -d'/' | cut -f2- -d'.' | rev)
mkdir $f
cd $f
module load  StdEnv/2020 gcc/9.3.0 sra-toolkit/2.10.8 htslib
source /project/6005588/shared/virtualenv_python_2.7/bin/activate
module load bwa samtools

tg=/project/6005588/shared/bin/trim_galore/TrimGalore-0.6.0/trim_galore
ivar_path=/lustre06/project/6005588/shared/bin/iVar/ivar-1.3/bin

amplicons_path=/lustre06/project/6005588/shared/References/covid-19/artic_ncov2019_primers/hybrid_amplicon_artic_v3_v4_v4.1/nCOV-2019.scheme.hybrid_v3_v4_v4.1.modCHRname.bed

ref=/project/6005588/shared/References/covid-19/NCBI/NC_045512.2/sequence.fasta
gff=/project/6005588/shared/References/covid-19/NCBI/NC_045512.2/GCF_009858895.2_ASM985889v3_genomic.gff

pe_or_se=$(fastq-dump -X 1 -Z --split-spot $file 2> /dev/null | wc -l | awk '$1!=8{print f,"single-end"}$1==8{print f,"paired-end"}' 2> /dev/null)

if [ $pe_or_se == "paired-end" ]; then 
fastq-dump -I --gzip --split-e $file

gunzip -c ${f}_1.fastq.gz | sed -E 's/(^[@+][ESD]RR[0-9]+\.[0-9]+)\.[12]/\1/' | gzip -c > ${f}_1_fixed.fastq.gz
gunzip -c ${f}_2.fastq.gz | sed -E 's/(^[@+][ESD]RR[0-9]+\.[0-9]+)\.[12]/\1/' | gzip -c > ${f}_2_fixed.fastq.gz

$tg --paired -q 20 -o $f.trimmedQ20 ${f}_1_fixed.fastq.gz ${f}_2_fixed.fastq.gz

bwa mem -o $f.bwa.sam $ref ${f}.trimmedQ20/${f}_1_fixed_val_1.fq.gz ${f}.trimmedQ20/${f}_2_fixed_val_2.fq.gz

samtools view -Sb $f.bwa.sam > $f.bwa.bam
samtools sort $f.bwa.bam > $f.bwa.sort.bam
samtools index $f.bwa.sort.bam
rm $f.bwa.bam $f.bwa.sam

${ivar_path}/ivar trim -i $f.bwa.sort.bam -b $amplicons_path -q 15 -m 30 -s 4 -p $f.bwa.amplicon_trim_smooth -e > $f.ivar.out 

samtools sort $f.bwa.amplicon_trim_smooth.bam > $f.bwa.amplicon_trim_smooth.sorted.bam
samtools view -F 2048 -bo $f.bwa.amplicon_trim_smooth.sorted.primaryOnly.bam $f.bwa.amplicon_trim_smooth.sorted.bam
samtools index $f.bwa.amplicon_trim_smooth.sorted.primaryOnly.bam
rm $f.bwa.amplicon_trim_smooth.sorted.bam

samtools mpileup -d 10000 -A -Q 0 $f.bwa.amplicon_trim_smooth.sorted.primaryOnly.bam | ${ivar_path}/ivar consensus -p $f.bwa.amplicon_trim_smooth.primaryOnly.consensus -q 20 -t 0.75 -m 20 > $f.samtools.mpileup.out


##########################
############# QC PART
##########################
module load java
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles

module load mugqic/qualimap/2.2.1

qualimap bamqc -nt 1 -outformat PDF -bam $f.bwa.sort.bam -outdir initialBam_$f --java-mem-size=10G

qualimap bamqc -nt 1 -outformat PDF -bam $f.bwa.amplicon_trim_smooth.sorted.primaryOnly.bam -outdir filteredBam_$f --java-mem-size=10G

samtools flagstat $f.bwa.amplicon_trim_smooth.sorted.primaryOnly.bam > $f.samtools_flagstat.out


(
#Library name
echo "Library_name" $f

#Number of initial reads
cat $f.trimmedQ20/${f}_1_fixed.fastq.gz_trimming_report.txt | grep "sequences processed in total" | awk '{print "input_R1",$1}'
cat $f.trimmedQ20/${f}_2_fixed.fastq.gz_trimming_report.txt | grep "sequences processed in total" | awk '{print "input_R2",$1}'

#Number of reads having adapters
cat $f.trimmedQ20/${f}_1_fixed.fastq.gz_trimming_report.txt | grep "Reads with adapters:" | awk '{print "tg_adapters_R1",$4}'
cat $f.trimmedQ20/${f}_2_fixed.fastq.gz_trimming_report.txt | grep "Reads with adapters:" | awk '{print "tg_adapters_R2",$4}' 

#Number of reads trimmed because of bad quality
cat $f.trimmedQ20/${f}_1_fixed.fastq.gz_trimming_report.txt | grep "Quality-trimmed:" | awk '{print "tg_badQual_R1",$2}' 
cat $f.trimmedQ20/${f}_2_fixed.fastq.gz_trimming_report.txt | grep "Quality-trimmed:" | awk '{print "tg_badQual_R2",$2}'

#Number of pairs removed because of one of the read shorter than the threshold (20)
cat $f.trimmedQ20/${f}_2_fixed.fastq.gz_trimming_report.txt | grep "the length cutoff " | awk '{print "tg_too_short",$19}'

#Test if the index exist for the bam after bwa and the final one
ls $f.bwa.sort.bam.bai*  | awk 'END{print NR==1 ? "alignment_index_exist YES" : "alignment_index_exist NO"}'
ls $f.bwa.amplicon_trim_smooth.sorted.primaryOnly.bam*  | awk 'END{print NR==1 ? "finalBAM_index_exist YES" : "finalBAM_index_exist NO"}'

#number of read and trimming 
cat initialBam_$f/genome_results.txt | grep "number of mapped reads" | awk '{print "mapped_reads_raw",$6}'
cat initialBam_$f/genome_results.txt | grep "number of supplementary alignments" | awk '{print "supp_alignement",$6}'
cat $f.ivar.out | grep "Trimmed primers from" | awk '{print "reads_trimmed_with_primers",substr($5,2,length($5)-2)}'
cat filteredBam_$f/genome_results.txt | grep "number of mapped reads" | awk '{print "primary_trimmed_mapped_reads",$6}'
cat filteredBam_$f/genome_results.txt | grep "first in pair" | awk '{print "mapped_reads_R1",$10}'
cat filteredBam_$f/genome_results.txt | grep "second in pair" | awk '{print "mapped_reads_R2",$10}'
cat filteredBam_$f/genome_results.txt | grep "number of duplicated reads" | awk '{print "duplicated_reads",$7}'
cat $f.samtools_flagstat.out | grep "properly paired"  | awk '{print "properly_paired",$1}'
cat filteredBam_$f/genome_results.txt | grep "mean insert size" | awk '{print "mean_insert_size",$5}'
cat filteredBam_$f/genome_results.txt | grep "median insert size" | awk '{print "median_insert_size",$5}'
cat filteredBam_$f/genome_results.txt | grep "number of A's" | awk '{print "number_of_A",$5}'
cat filteredBam_$f/genome_results.txt | grep "number of C's" | awk '{print "number_of_C",$5}'
cat filteredBam_$f/genome_results.txt | grep "number of G's" | awk '{print "number_of_G",$5}'
cat filteredBam_$f/genome_results.txt | grep "number of T's" | awk '{print "number_of_T",$5}'
cat filteredBam_$f/genome_results.txt | grep "mapped reads with deletion percentage" | awk '{print "percent_mapped_reads_with_DEL",$7}'
cat filteredBam_$f/genome_results.txt | grep "mapped reads with insertion percentage" | awk '{print "percent_mapped_reads_with_INS",$7}'
cat filteredBam_$f/genome_results.txt | grep "mean coverageData" | awk '{print "meann_coverage",substr($4,1,length($4)-1)}'
cat filteredBam_$f/genome_results.txt | grep "reference with a coverageData >= 10X" | awk '{print "percent_genome_covered_over_10X",$4}'

) | sed 's/[,%]//g'  | tr ' ' '\t' > $f.qc.tsv

cp $f.bwa.amplicon_trim_smooth.sorted.primaryOnly.bam ..
cp $f.bwa.amplicon_trim_smooth.sorted.primaryOnly.bam.bai .. 
cp $f.qc.tsv ..
cp $f.bwa.amplicon_trim_smooth.primaryOnly.consensus ..
rm *fastq.gz
rm SRR13268199.trimmedQ20/*fq.gz
cd ..
tar -czvf $f.tar.gz $f
rm -r $f
fi
echo $f is $pe_or_se
