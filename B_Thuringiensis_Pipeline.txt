#!/bin/bash

# Download the paired end reads for B thuringiensis

echo -e "Initializing Bioinformatics Pipeline...\n...\n...\n..."

echo "Retrieving fastqc file with wget..."
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/001/SRR2093871/SRR2093871_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/001/SRR2093871/SRR2093871_2.fastq.gz
echo "Files retrieved"

# Run FastQC
echo "Analyzing fastqc file..."
fastqc SRR2093871_1.fastq.gz  
fastqc SRR2093871_2.fastq.gz 
echo "Analysis complete"

# fastqc should indicate no problems - move on to assembly
cd ../
mkdir 01_spades
cd 01_spades

# Run spades with four threads and the pair-end reads (use half your total threads available)

echo "\nAssembling seqeunce reads..."
spades.py \
-k 21,33,55,77,99,127 \
-t 4 \
--pe1-1 ../00_rawdata/SRR2093871_1.fastq.gz  \
--pe1-2 ../00_rawdata/SRR2093871_2.fastq.gz  \
--careful \
-o assembly
echo "Assembly completed"

# Now let's use quast to test our scaffolds

./quast.py -t 4 ~/BIFS619/Bthuringiensis/01_spades/assembly/scaffolds.fasta

echo "Viewing Results..."
# Examine the results:
cd quast_results/latest/
ff report.html

#########################################################################################

# Moving on to Gene Prediction and Annotation

echo -e "\nGene Annotation Initializing...\n...\n"
cd ~/BIFS619/Bthuringiensis/02_annotations

cp ~/BIFS619/Bthuringiensis/01_spades/assembly/scaffolds.fasta .

echo "Running Prodigal..."
prodigal -i scaffolds.fasta \
-o prodigal.out -s prodigal_table.out -f gff -d prodigal_geneseqs.out \
-a prodigal_translated.out

echo "Running Glimmer..."
tigr-glimmer long-orfs -n -t 1.15 scaffolds.fasta glimmer.longorf

# creating the training set for long ORFs
tigr-glimmer extract -t scaffolds.fasta glimmer.longorfs > glimmer.train

# fitting an HMM model to our training data
tigr-glimmer build-icm -r glimmer.icm < glimmer.train

# predicting genes using the HMM model that was just created
tigr-glimmer glimmer3 -o50 -g110 -t30 scaffolds.fasta glimmer.icm glimmer

echo "Running Barrnap..."
barrnap --threads 4 --kingdom bac scaffolds.fasta > rrna.gff

conda init bash

conda activate prokka_env

prokka scaffolds.fasta

#########################################################################################

# Gene Annotation Output Parsing

echo -e "\nInitializing Annotation Output Parsing...\n...\n...\n..."

# Adding Barrnaap output to final.gff file
sed 1d rrna.gff | sed 's/_length_.*_cov_.*//g' > rrna_id.tmp

sed 1d rrna.gff | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > rrna_data.tmp

touch final.gff

echo -e "seqID\tsource_app\tfeature_type\tfeat_start\tfeat_stop\tscore\tstrand\tphase\tattributes" >> final.gff

paste rrna_id.tmp rrna_data.tmp >> final.gff

# removing temporary files
rm *.tmp

# Adding Glimmer results to final .gff

sed 1,22d glimmer.detail > glimmer.1.tmp

# breaking up glimmer files into separate files for each scaffold
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.temp
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < glimmer.1.tmp


# filtering through glimmer outfiles the best gene hits
for file in ls *.temp
do
        echo "Now Processing $file"
        # Find lines that begin with numbers
        awk --posix '{ if ($1 ~ /^[0-9]/) print $1}' $file > temp_holding.tmp
                # if lines found, truncate $file to ID and add entries to gff
                if [ ! -s temp_holding.tmp ];
                        then
                                echo "No predicted genes found";
                        else
                                echo "Predicted Genes being processed";
                                prefix=$(echo $file | sed 's/_length.*//')
                                while read line
                                do
                                        # $line contains the ID number we should parse out of $file
                                        # need to truncate $file
                                        ID=$(echo $file | sed 's/_length.*//')                                  
                                        grep "^$line" $file | awk '{print $0}' | awk -v ID="$ID" \
                                        '{print ID, "glimer", "CDS", $4, $5, $9, $2, ".", "."}' >> final.gff
                                done < temp_holding.tmp
                fi
done

# Parsing Prodigal's output and adding to final.gff

grep -v "#" prodigal.out > prodigal.tmp.gff

# Generate ID file
sed 's/_length_.*//g' prodigal.tmp.gff > prodigal.id.tmp

# Generate "everything else" file
awk '{print $2, $3, $4, $5, $6, $7, $8, $9}' prodigal.tmp.gff > prodigal.data.tmp

# paste the stuff together
paste prodigal.id.tmp prodigal.data.tmp > prodigal.gff

# Add to the final.gff table
cat prodigal.gff >> final.gff

# clean up the format a bit
sed 's/_v2.6.3//g' final.gff > final.gff.tmp
cat final.gff.tmp | tr ' ' '\t' > final.gff.tmp2
rm final.gff
mv final.gff.tmp2 final.gff

# remove tmp files
rm *tmp*

#########################################################################################

# RNA Seq

echo -e "\nInitializing RNASeq and the Transcriptome Analysis...\n...\n..."
cd ~/BIFS619/Bthuringiensis/03_RNASEQ
echo "Retrieving RNASeq files in the .gz format..."
wget https://culture-bioinformatics.org/UMGC/BIFS619/SRR2093871_1.fastq.gz

wget https://culture-bioinformatics.org/UMGC/BIFS619/SRR2093871_2.fastq.gz

echo "Decompressing .gz files..."
gunzip *.gz

echo "Retrieving our reference genome FASTA file..."
wget https://culture-bioinformatics.org/UMGC/BIFS619/Bthuringiensis.fasta

echo "Mapping the first file..."
minimap2 -ax sr Bthuringiensis.fasta SRR2093871_1.fastq > Bt1.sam
echo "Map 1 complete"

echo "Mapping the second file"
minimap2 -ax sr Bthuringiensis.fasta SRR2093871_2.fastq > Bt2.sam

echo "Converting sam files to bam files"
samtools view -b Bt1.sam > Bt1.bam
samtools view -b Bt2.sam > Bt2.bam

echo "Sorting the bam files"
samtools sort Bt1.bam -o Bt1.sorted.bam
samtools sort Bt2.bam -o Bt2.sorted.bam

echo "Retrieving Gene Map..."
wget https://culture-bioinformatics.org/UMGC/BIFS619/Bthuringiensis_gene_map.txt

echo "Running bam-readcount 1..." 
bam-readcount -w 0 -q 0 -b 0 -f Bthuringiensis.fasta Bt1.sorted.bam > counts.1.txt

echo "Running bam-readcount 1..."
bam-readcount -w 0 -q 0 -b 0 -f Bthuringiensis.fasta Bt2.sorted.bam > counts.2.txt

echo "Running counts1..."
awk '{print $2, $4}' counts.1.txt > simple.counts.1.txt

echo "Running counts2..."
awk '{print $2, $4}' counts.2.txt > simple.counts.2.txt


touch genecounts1.txt
while read id start stop mid gene;
do
	count=$(grep -w $mid simple.counts.1.txt | awk '{print $2}')
if [[ "$count" == "" ]] ; then 
  count=0
fi
	echo -e "$count\t $gene" >> genecounts1.txt

done < Bthuringiensis_gene_map.txt

touch genecounts2.txt
while read id start stop mid gene;
do
        count=$(grep -w $mid simple.counts.2.txt | awk '{print $2}')
if [[ "$count" == "" ]] ; then
  count=0
fi
        echo -e "$count\t $gene" >> genecounts2.txt

done < Bthuringiensis_gene_map.txt

echo -e "\n\n\nScript complete!"








