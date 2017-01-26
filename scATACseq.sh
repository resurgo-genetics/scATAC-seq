#PIPELINE FOR scATAC-seq ALIGNMENT AND ANALYSIS
#have to add something somewhere to remove the header in the first place

	#alignment
	
		#build bwa index file
		bwa index -p hg19bwaidx -a bwtsw /home/cnb3001/GRCh37.primary_assembly.genome.fa.gz
		#align sequences
		bwa aln -t 4 hg19bwaidx SRR1947692_1.fastq > SRR1947692_1.bwa
		bwa aln -t 4 hg19bwaidx SRR1947692_2.fastq > SRR1947692_2.bwa
		bwa sampe hg19bwaidx SRR1947692_1.bwa SRR1947692_2.bwa SRR1947692_1.fastq SRR1947692_2.fastq > SRR1947692.sam
		
	#processing alignment
		#convert to bam
		samtools view -Shb SRR1947692.sam > SRR1947692.bam

		#remove poorly mapped reads and sort
		samtools view -bhq 10 SRR1947692.bam > SRR1947692_MQ.bam
		samtools sort SRR1947692_MQ.bam -o SRR1947692_MQ.sorted.bam > SRR1947692_MQ.sorted.bam
		
		#remove duplicates
		java -jar /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/kwe2001/programs/picard-tools-1.85/MarkDuplicates.jar INPUT=SRR1947692_MQ.sorted.bam OUTPUT=SRR1947692_MQ_dedup.sorted.bam METRICS_FILE=SRR1947692_dups.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true

		samtools view -h SRR1947692_MQ_dedup.sorted.bam > SRR1947692_MQ_dedup.sorted.sam
		
		grep -v "chrM" SRR1947692_MQ_dedup.sorted.sam > SRR1947692_MQ_dedup_noM.sorted.sam

	
	#call peaks on aggregate signal using homer, store in bed file
		mkdir tagdir
		#may have to convert back to sam first, check this
		makeTagDirectory tagdir SRR1947693_MQ_dedup_noM.sorted.sam
		findPeaks tagdir -o SRR1947692_peaks.out
		grep -v "#" SRR1947692_peaks.out | cut -f2-5 > SRR1947692_peaks.bed

	#deriving count matrix for each cell; this is a custom python script
		#run first part of script
		python /Users/cass/Documents/LeslieLab/Git_Repo/scATAC-seq/scATAC_peakmatrix_cl1.py -s SRR1947692_MQ_dedup_noM.sorted.sam -b GSM1647122_GM12878vsHEK.indexconversion.txt -o test.sam -r SRR1947692
		#add header back, convert to bam, and index output
		#when using older version of samtools, may need to add a -S option for view command
		cat header.txt test.sam > test_header.sam
		samtools view -b test_header.sam > test_header.bam
		samtools index test_header.bam
		
		#run second part of script
		python /Users/cass/Documents/LeslieLab/Code/scATAC_peakmatrix_cl2.py -b test_header.bam -p SRR1947692_peaks_new.bed -n 38995
	#clustering in R