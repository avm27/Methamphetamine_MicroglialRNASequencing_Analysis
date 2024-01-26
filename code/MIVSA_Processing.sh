#!/bin/bash

python="/home/avm27/anaconda3/envs/RNA_Seq/bin/python"
stringtiescript="/home/avm27/anaconda3/envs/RNA_Seq/bin/prepDE.py"
experimentname="MIVSA"
experimentnamefull="Methamphetamine Intravenous Self Administration"

echo "Final Code: MIVSA Experiment Pegasus Pre-Processing"

echo "Loading Environment"

cd /home/avm27/Documents/Raw_Sequencing_Data/MIVSA

for histName in $(find /home/avm27/Documents/Raw_Sequencing_Data/MIVSA -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
do 

    echo "Merging R1" ${histName}

    cat ${histName}_L00*_R1_001.fastq.gz > ${histName}_R1.fastq.gz
   
    echo "Merging R2" ${histName}

    cat ${histName}_L00*_R2_001.fastq.gz > ${histName}_R2.fastq.gz

done
wait

echo "Trimming Samples"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/MIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    trim_galore --cores 6 --dont_gzip --paired ${samplename}_R1.fastq.gz ${samplename}_R2.fastq.gz

wait
done
wait

echo "Trimming Completed"

echo "Begin Aligning to the mm10 Genome"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/MIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    STAR --runThreadN 8 \
    --runMode alignReads  \
    --genomeDir /home/avm27/genome/mm10/ncbi_STAR \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD XS \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNoverLmax 0.3 \
    --outFileNamePrefix ${samplename} \
    --quantMode TranscriptomeSAM GeneCounts \
    --readFilesIn ${samplename}_R1_val_1.fq ${samplename}_R2_val_2.fq \
    --outTmpDir ${samplename}

wait

gzip ${samplename}_R1_val_1.fq
gzip ${samplename}_R2_val_2.fq

wait
done
wait


echo "Obtain unique mapped reads"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/MIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    samtools view ${samplename}Aligned.sortedByCoord.out.bam | grep -w 'NH:i:1' | samtools view -bt /home/avm27/genome/mm10/mm10.fa.fai > ${samplename}.bam

done
wait


echo "Index Samples"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/MIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    samtools index ${samplename}.bam

wait
done
wait

echo "Create GTF Files"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/MIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    stringtie -p 12 -G /home/avm27/genome/mm10/mm10.ncbiRefSeq.gtf -e -o ${samplename}.gtf -A ${samplename}_FPKM.txt -l ${samplename} ${samplename}.bam

wait
done


## Generate Sample Info List for Stringtie. Name the file "sample_lst.txt" and save to the working directory. Below is the code ran for this experiment.

echo "Generating sample info list for StringTie"

echo "SAL01 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_SAL01.gtf" > sample_lst_$experimentname.txt
echo "SAL02 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_SAL02.gtf" >> sample_lst_$experimentname.txt
echo "SAL03 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_SAL03.gtf" >> sample_lst_$experimentname.txt
echo "SAL04 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_SAL04.gtf" >> sample_lst_$experimentname.txt
echo "SAL05 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_SAL05.gtf" >> sample_lst_$experimentname.txt
echo "SAL06 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_SAL06.gtf" >> sample_lst_$experimentname.txt
echo "SAL07 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_SAL07.gtf" >> sample_lst_$experimentname.txt
echo "SAL08 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_SAL08.gtf" >> sample_lst_$experimentname.txt
echo "MN01 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_MN01.gtf" >> sample_lst_$experimentname.txt
echo "MN02 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_MN02.gtf" >> sample_lst_$experimentname.txt
echo "MN03 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_MN03.gtf" >> sample_lst_$experimentname.txt
echo "MN04 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_MN04.gtf" >> sample_lst_$experimentname.txt
echo "MN05 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_MN05.gtf" >> sample_lst_$experimentname.txt
echo "MN06 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_MN06.gtf" >> sample_lst_$experimentname.txt
echo "MN07 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_MN07.gtf" >> sample_lst_$experimentname.txt
echo "MN08 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_MN08.gtf" >> sample_lst_$experimentname.txt
echo "CRV01 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_CRV01.gtf" >> sample_lst_$experimentname.txt
echo "CRV02 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_CRV02.gtf" >> sample_lst_$experimentname.txt
echo "CRV03 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_CRV03.gtf" >> sample_lst_$experimentname.txt
echo "CRV04 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_CRV04.gtf" >> sample_lst_$experimentname.txt
echo "CRV05 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_CRV05.gtf" >> sample_lst_$experimentname.txt
echo "CRV06 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_CRV06.gtf" >> sample_lst_$experimentname.txt
echo "CRV07 /home/avm27/Documents/Raw_Sequencing_Data/MIVSA/MIVSA_CRV07.gtf" >> sample_lst_$experimentname.txt

echo "Completed Same Info File"

## Generate Gene Count matrix for DESEQ2. Name the file "gene_count_matrix_MIVSA.csv" and save to the working directory.

echo "Generating Gene Count Matrix File"

$python $stringtiescript -i sample_lst.txt -g gene_count_matrix_MIVSA.csv

echo $experimentnamefull "RNA-Sequencing Pre-Processing Analysis Complete"

