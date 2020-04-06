#start with gzipped fastqs in a folder named "fastq" and a .txt file named "filePairs" which contains R1 and R2 reads separated by ";" in the fastq folder (no.gz)
source /etc/profile.d/apps.sh

echo "Program Started: $(date)" > timelog.txt

parallel gunzip ::: ./fastq/*.gz

echo "Files unzipped: $(date)" >> timelog.txt

cd fastq/

for FILEPAIR in $(cat filepairs.txt)
do
	F1=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $1}')
	F2=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $2}')
	cutadapt -a CTGTCTCTTATACACATCT \
		-A AGATGTGTATAAGAGACAG \
		--minimum-length=25 \
		-j 16 \
		-o "trimmed_${F1}"  -p "trimmed_${F2}" \
		${F1} ${F2} &
done
wait

cd ../

echo "Adapters Trimmed: $(date)" >> timelog.txt

mkdir trimmedFastq
mv ./fastq/filepairs.txt ./trimmedFastq

for FILE in ./fastq/trimmed*
	do
	   {
		mv $FILE ./trimmedFastq
	   } &
	done
wait


mkdir BAM

cd trimmedFastq/

for FILEPAIR in $(cat filepairs.txt)

do
	F1=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $1}')
	F2=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $2}')
	echo "beginning alignment to galGal6 of file $FILEPAIR"
	bowtie2 --local --very-sensitive-local \
		--no-unal --no-mixed --no-discordant \
		--threads 16 \
		-x /home/ash274/workdir/genome/bowtie2-GRCg6a/GRCg6a \
		-X 2000 \
		-1 "trimmed_${F1}" \
		-2 "trimmed_${F2}" \
		-S ${F1/_R1.fastq/toGalGal6.SAM}
	echo "complete"
done
wait

cd ../

echo "Reads Aligned: $(date)" >> timelog.txt


for FILE in ./trimmedFastq/*.SAM
do
	samtools view -bS -@ 16 $FILE |
    samtools sort -@ 16 > "${FILE/.SAM/.BAM}" &&
    samtools index -@ 16 "${FILE/.SAM/.BAM}"
done
wait

for FILE in ./trimmedFastq/*.BAM
	do
	{
		mv $FILE ./BAM
	} &
	done
wait

for FILE in ./trimmedFastq/*.bai
	do
	{
		mv $FILE ./BAM
	} &
	done
wait

for file in ./BAM/*.BAM

do
	java -jar $PICARD MarkDuplicates \
			I="${file}" \
			O="${file/.BAM/_dupsMarked.BAM}" \
			M="${file/.BAM/_dupsMarkedStats.txt}"
done
wait

echo "Duplicates marked: $(date)" >> timelog.txt

for FILE in ./BAM/*_dupsMarked.BAM
do
	samtools view -@ 16 -F 1804 -f 2 -b $FILE |
	samtools sort -@ 16 > "${FILE/_dupsMarked.BAM/_nodups.BAM}" &&
	samtools index "${FILE/_dupsMarked.BAM/_nodups.BAM}"
done
wait

echo "Duplicates removed: $(date)" >> timelog.txt

for FILE in ./BAM/*_nodups.BAM
do
	bamCoverage --bam "$FILE" \
			--outFileName "$FILE".bw \
			--outFileFormat bigwig \
			--binSize 5 \
			--numberOfProcessors 15 \
			--normalizeUsing RPGC \
			--effectiveGenomeSize 1218492533 \
			--extendReads		
done
wait

echo "BigWigs made: $(date)" >> timelog.txt

mkdir BW

for FILE in ./BAM/*.bw
	do
	{
		mv $FILE ./BW
	} &
	done
wait

for FILE in ./BAM/*_nodups.BAM

do 
	macs2 callpeak -t "$FILE" \
			-n "$FILE" \
			-f BAMPE \
			-g 1218492533 \
			-q 0.05 \
			--call-summits \
			--nomodel --shift 37 \
			--ext 73 -B --SPMR
done