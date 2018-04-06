#############################################################
### F1 Bioinformatics 
### Practical
### Part: Metabarcoding
### Lecturer: PD Dr. Alexander Keller
#############################################################
# This is material of a guided lecture. 
# Individual commands are explained while progressing through the script
# outputs are interpreted together on occurrence
# as well as intermediate files checked with "less" or "head"
#############################################################


#############################################################
# create a new directory to work in
mkdir -p F1_MBC
cd F1_MBC

#############################################################
#### Obtaining data
# Getting the Sequence Data Information
wget http://www.ebi.ac.uk/ena/data/warehouse/filereport\?accession\=PRJEB8640\&result\=read_run\&fields\=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy\&download\=txt -O reads.tsv

# we only want 20 of these
head reads.tsv > reads_subset.tsv
tail reads.tsv >> reads_subset.tsv
less -S reads_subset.tsv
#press q to exit the text reader

# Download the Raw Data
for i in $(cut -f13 reads_subset.tsv | grep fastq.gz | perl -pe 's/;/\n/')
do
    wget $i
done

### Stop here for the introduction talk, the download takes a couple of minutes

#############################################################
### Setting variables
data=$(pwd)
u9=/storage/full-share/mbd/usearch9                     # USEARCH 9.0 binary
u8=/storage/full-share/mbd/usearch8                     # USEARCH 8.0 binary
f=/storage/full-share/mbd/fastq-join/fastq-join         # FASTQ-JOIN binary
p=/storage/full-share/mbd/python_scripts                # Folder of python scripts
s=/storage/full-share/mbd/seqfilter/bin/SeqFilter       # SeqFilter binary

#############################################################
### Setting up the analysis
mkdir -p $data/raw
mkdir -p $data/joined

gunzip *.gz
ls -l

#store file suffixes as variables for forward and reverse reads

RF='_R1_001.fastq';
RR='_R2_001.fastq';

# get sample names
ls  $data/*$RF | sed "s/^.*\/\([a-zA-Z0-9_.-]*\)$/\1/g" | sed "s/$RF//" > samples.txt
head samples.txt



#############################################################
### Data processing

# Read Joining and Quality Filtering
for file in `cat samples.txt` ;
do
 	echo "Processing >>> $file <<<";
	echo "..join ends";

 	# joining forward and reverse reads
	$f $data/$file$RF $data/$file$RR -o $data/joined/$file.%.fq

	# keep R1 files for those that do not join, perhaps they are long enough and of good quality
	cat $data/joined/$file.join.fq $data/joined/$file.un1.fq > $data/joined.$file.fastq

	# move original files
	mv $data/$file$RF $data/raw/; mv $data/$file$RR $data/raw/

	# filter reads with high expected error rate, that are too short or have Ns
 	echo "..filter";
 	$u9 -fastq_filter $data/joined.$file.fastq -fastq_maxee 1 -fastq_minlen 200 -fastq_maxns 1 -fastaout filter.$file.fasta -threads 4

	# rename sequence names to match their original sample
 	echo "..parse";
	python $p/fasta_number.py filter.$file.fasta $file. > parsed1.$file.fasta
	cat parsed1.$file.fasta | sed "s/_L001//g" | sed "s/\./_/g " > parsed2.$file.fasta

done

### a good time for questions now!

#############################################################
### Taxonomic classification

### Direct Hits
# Get reference database
wget https://github.com/molbiodiv/meta-barcoding-dual-indexing/raw/master/data/viridiplantae_bavaria_2015.fa
head viridiplantae_bavaria_2015.fa

# Matching
# combine files of all samples to a single file to be searched
cat parsed2.* > all.fasta

 
# clean up temporary files
rm -r raw/ joined/ parsed* joined.* filter.*

# convert to a barcoding format readable by usearch
cat all.fasta | sed -e "s/^>\([a-zA-Z0-9-]*\)_\(.*$\)$/>\1_\2;barcodelabel=\1/" > all.bc.fasta

# find best direct hits for all filtered sequences with more than 97% identity
$u9 -usearch_global all.bc.fasta -db viridiplantae_bavaria_2015.fa -id 0.97 -uc output_BV3.uc -fastapairs output_BV3.fasta  -strand plus -threads 4

### Hierarchic Classification
# Get reference database
wget https://github.com/iimog/meta-barcoding-dual-indexing/raw/v1.1/training/utax/utax_trained.tar.gz

tar xzvf utax_trained.tar.gz
rm utax_trained.tar.gz

$u8 -makeudb_usearch utax_trained/viridiplantae_all_2014.utax.fa -output utax_trained/viridiplantae_all_2014.utax.udb

# get names and sequences of those without hits, classify them herarchical with utax afterwards

grep "^N[[:space:]]" output_BV3.uc | cut -f 9 > output_BV3.nohit
$s all.bc.fasta --ids output_BV3.nohit --out all.bc.BV3.nohit.fasta
$u8 -utax all.bc.BV3.nohit.fasta -db utax_trained/viridiplantae_all_2014.utax.udb -utax_rawscore -tt utax_trained/viridiplantae_all_2014.utax.tax -utaxout all.bc.BV3.nohit.utax

### another good time for questions now!

less -S all.bc.BV3.nohit.utax

# create a pseudo.uc file of the hierarchical classification and filter low quality assignments
perl -ne '
($id, $tax, $sign)=split(/\t/);
@tmp=split(/,/, $tax);
@tax=();
foreach $t (@tmp){
    $t=~s/__/:/g;
    $t=~s/ +/_/g;
    $t=~s/\(([\d.]+)\)//;
    if($1 < 27){last;}
    push @tax, $t;
    $t=~/_(\d+)$/;
    $taxid=$1;
}
next if(@tax == 0);
print "H".("\t"x8)."$id\tt$taxid;tax=".join(",", @tax).";\n"
' all.bc.BV3.nohit.utax | sed "s/,s:.*$//g" > all.bc.BV3.nohit.pseudo.uc

less -S all.bc.BV3.nohit.pseudo.uc

#############################################################
### Preparing Data for R

#combine direct hit .uc and the hierarchically classified pseudo.uc
cat output_BV3.uc all.bc.BV3.nohit.pseudo.uc > combined.uc

#convert output format to TAX/OTU-Table (community matrix including taxonomic lineages)
python $p/uc2otutab.py  combined.uc > combined.txt

less -S combined.txt

#split TAX-OTU-Table to OTU and TAX Table respectively
cat combined.txt | sed "s/;tax=.*;//g;s/:/_/g" > combined.otu
cat combined.txt | cut -f 1 | sed "s/;tax=/,/g;s/:/_/g" | sed "s/;//" | sed "s/OTUId/,Kingdom,Phylum,Class,Order,Family,Genus,Species/" > combined.tax

head combined.otu
head combined.tax
tail combined.tax

cat reads_subset.tsv | cut -f8,13 | sed -e "s/ftp.sra.ebi.ac.uk.*\///" -e "s/_.*//" -e "s/\(.*\)\t\(.*\)/\2\t\1/" > MapFile.txt 
sed  -i "s/scientific/#SampleID\tBeeSpecies/" MapFile.txt

less MapFile.txt



# now proceed in R with the F1_MBC.R
