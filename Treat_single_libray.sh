file=$1
f=$(echo $file | rev | cut -f1 -d'/' | cut -f2- -d'.' | rev)
module load  StdEnv/2020 gcc/9.3.0 sra-toolkit/2.10.8 htslib
source /project/6005588/shared/virtualenv_python_2.7/bin/activate
module load bwa samtools

tg=/project/6005588/shared/bin/trim_galore/TrimGalore-0.6.0/trim_galore
ivar_path=/lustre06/project/6005588/shared/bin/iVar/ivar-1.3/bin

amplicons_path=/lustre06/project/6005588/shared/References/covid-19/artic_ncov2019_primers/hybrid_amplicon_artic_v3_v4_v4.1/nCOV-2019.scheme.hybrid_v3_v4_v4.1.modCHRname.bed

ref=/project/6005588/shared/References/covid-19/NCBI/NC_045512.2/sequence.fasta
gff=/project/6005588/shared/References/covid-19/NCBI/NC_045512.2/GCF_009858895.2_ASM985889v3_genomic.gff

pe_or_se=$(fastq-dump -X 1 -Z --split-spot $file 2> /dev/null | wc -l | awk '$1!=8{print f,"single-end"}$1==8{print f,"paired-end"}' 2> /dev/null)

if [ "$pe_or_se" == "paired-end" ]; then 
fastq-dump -I --gzip --split-files $file.sra

gunzip -c ${f}_1.fastq.gz | sed -E 's/(^[@+][ES]RR[0-9]+\.sra\.[0-9]+)\.[12]/\1/' | gzip -c > ${f}_1_fixed.fastq.gz
gunzip -c ${f}_2.fastq.gz | sed -E 's/(^[@+][ES]RR[0-9]+\.sra\.[0-9]+)\.[12]/\1/' | gzip -c > ${f}_2_fixed.fastq.gz

$tg --paired -q 20 -o $f.trimmedQ20 ${f}_1_fixed.fastq.gz ${f}_2_fixed.fastq.gz

bwa mem -o $f.bwa.sam $ref ${f}.trimmedQ20/${f}_1_fixed_val_1.fq.gz ${f}.trimmedQ20/${f}_2_fixed_val_2.fq.gz

samtools view -Sb $f.bwa.sam > $f.bwa.bam
samtools sort $f.bwa.bam > $f.bwa.sort.bam
samtools index $f.bwa.sort.bam
rm $f.bwa.bam $f.bwa.sam

${ivar_path}/ivar trim -i $f.bwa.sort.bam -b $amplicons_path -q 15 -m 30 -s 4 -p $f.bwa.sort.amplicon_trim_smooth -e 

samtools sort $f.bwa.sort.amplicon_trim_smooth.bam > $f.bwa.amplicon_trim_smooth.sorted.bam
samtools index $f.bwa.amplicon_trim_smooth.sorted.bam

samtools mpileup -d 10000 -A -Q 0 $f.bwa.amplicon_trim_smooth.sorted.bam | ${ivar_path}/ivar consensus -p $f.bwa.amplicon_trim_smooth.consensus -q 20 -t 0.75 -m 20

#QC PART
nb_r1=$(zcat ${f}_1.fastq.gz | awk 'END{print NR/2}')
nb_r2=$(zcat ${f}_2.fastq.gz | awk 'END{print NR/2}')
libtype=$(grep $f logs/log.$mod.$step | egrep "is paired-end|is single" | cut -d' ' -f3)
nb_r1_ad=$(cat $f.trimmedQ20/${f}_1_fixed.fastq.gz_trimming_report.txt | grep "Reads with adapters:" | awk '{print $4}' | sed 's/,//g' )
nb_r2_ad=$(cat $f.trimmedQ20/${f}_2_fixed.fastq.gz_trimming_report.txt | grep "Reads with adapters:" | awk '{print $4}' | sed 's/,//g' )
nb_bp1_qual=$(cat $f.trimmedQ20/${f}_1_fixed.fastq.gz_trimming_report.txt | grep "Quality-trimmed:" | awk '{print $2}' | sed 's/,//g' )
nb_bp2_qual=$(cat $f.trimmedQ20/${f}_2_fixed.fastq.gz_trimming_report.txt | grep "Quality-trimmed:" | awk '{print $2}' | sed 's/,//g' )
nb_r1_tg=$(cat $f.trimmedQ20/${f}_1_fixed.fastq.gz_trimming_report.txt | grep "sequences processed in total" | awk '{print $1}')
nb_r2_tg=$(cat $f.trimmedQ20/${f}_2_fixed.fastq.gz_trimming_report.txt | grep "sequences processed in total" | awk '{print $1}')
nb_r_tg=$(cat $f.trimmedQ20/${f}_2_fixed.fastq.gz_trimming_report.txt | grep "Number of sequence pairs removed because at least one read was shorter than the length cutoff " | awk '{print $19}')
#nb_read_Sam=$(cat $f.bwa.sam | cut -c-100 | grep -v ^@ | wc -l)
nb_read_Bam=$(samtools view $f.bwa.sort.bam | cut -f 1 | sort -u | wc -l)
bai_exist=$(ls $f.bwa.sort.bam.bai*  | awk 'END{print NR==1 ? "YES" : "NO"}')
last_bam_lines=$(samtools view $f.bwa.amplicon_trim_smooth.sorted.bam  | cut -f 1 | sort -u | wc -l)
last_bai_exist=$(ls $f.bwa.amplicon_trim_smooth.sorted.bam.bai*  | awk 'END{print NR==1 ? "YES" : "NO"}')
#insert size
#properly paired (samtools flag stat)
#mean coverage (average read per postion)
#breath of coverage (nb of position under 100)
#breath of coverage (nb of position under 400)
echo $f $step $mod $nb_r1 $nb_r2 $libtype $nb_r1_ad $nb_r2_ad $nb_bp1_qual $nb_bp2_qual $nb_r1_tg $nb_r2_tg $nb_r_tg $nb_read_Bam $bai_exist $last_bam_lines $last_bai_exist | tr ' ' '\t' > $f.stat.tsv
rm $f.bwa.bam.bai $f.bwa.bam $f.bwa.sam $f.bwa.sort.bam $f.last_bam_lines $f.bwa.sort.bam.bai $f.bwa.amplicon_trim_smooth.bam ${f}_1_fixed.fastq.gz ${f}_2_fixed.fastq.gz ${f}_1.fastq.gz ${f}_2.fastq.gz ${f}.bwa.sort.amplicon_trim_smooth.bam
rm -r $f.trimmedQ20
fi
