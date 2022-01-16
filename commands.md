# Here are the commands used to replicate the benchmark analysis in Jiang et al. 2021

### Do note that there are several necessary modifications made because of errors encountered throughout the analysis.

### Tool versions
* VISOR (v1.1)
* NanoVar (v1.3.8)
* Sniffles (v1.0.12)
* Truvari (v3.0.1)

### Downloading data
```
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd162.GRCh37.variant_call.vcf.gz
gzip -d nstd162.GRCh37.variant_call.vcf.gz

wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd137.GRCh37.variant_call.vcf.gz
gzip -d nstd137.GRCh37.variant_call.vcf.gz

curl -s ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > human_hs37d5.fasta.gz
gunzip human_hs37d5.fasta.gz
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' human_hs37d5.fasta
```

### Preparation of SVs for integration (Modified)
```
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep DEL | awk -F '\t' '{print $1"\t"$2}' > del_col1.txt
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep DEL | awk -F '\t' '{print $8}' | awk -F ';END=' '{print $2}' | awk -F ';' '{print $1"\tdeletion\tNone\t0"}' > del_col2.txt
paste del_col1.txt del_col2.txt | uniq > del.bed

grep EXPERIMENT=2 nstd137.GRCh37.variant_call.vcf | grep INS |awk -F '\t' '{print $1"\t"$2-1"\t"$2}' > ins_col1.txt
grep EXPERIMENT=2 nstd137.GRCh37.variant_call.vcf | grep INS |awk -F '\t' '{print $8}' | awk -F ';SEQ=' '{print "insertion\t" $2 "\t0"}' >ins_col2.txt
paste ins_col1.txt ins_col2.txt | uniq > ins.bed  # Modified (The original output file was saved as "del.bed", which does not make sense.)

grep CHM1 nstd162.GRCh37.variant_call.vcf | grep DUP | awk -F '\t' '{print $1"\t"$2}' > dup_col1.txt
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep DUP | awk -F '\t' '{print $8}' | awk -F ';END=' '{print $2}' | awk -F ';' '{print $1"\ttandem duplication\t2\t0"}' > dup_col2.txt
paste dup_col1.txt dup_col2.txt | uniq > dup.bed

grep CHM1 nstd162.GRCh37.variant_call.vcf | grep INV | awk -F '\t' '{print $1"\t"$2}' > inv_col1.txt
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep INV | awk -F '\t' '{print $8}' | awk -F ';END=' '{print $2}' | awk -F ';' '{print $1"\tinversion\tNone\t0"}' > inv_col2.txt
paste inv_col1.txt inv_col2.txt | uniq > inv.bed

cat del.bed ins.bed dup.bed inv.bed | sort -k 1,1 -k 2,2n > chm1.bed
```

### Integrate SVs using VISOR (Modified)
```
VISOR HACk -g ../resource/human_hs37d5.fasta -b ../resource/chm1.bed -o chm1  # Modified (The original command line raised an error as the "SHORtS.LASeR.bed file" should not be used at this step.)
```

### Simulate long reads (Modified)
#### The "SHORtS.LASeR.bed" file provided was modified to replace the 0s of start coordinates with 1s (Please refer to https://github.com/davidebolo1993/VISOR/issues/18).
```
# 3X Coverage
VISOR LASeR -g ../resource/human_hs37d5.fasta -s chm1/ -b ../resource/SV_evaluation/SHORtS.LASeR.edit.bed -o 3x_20k_90 --coverage 3 --identity_min 90 --length_mean 20000 --read_type pacbio --threads 40
#### 5X Coverage
VISOR LASeR -g ../resource/human_hs37d5.fasta -s chm1/ -b ../resource/SV_evaluation/SHORtS.LASeR.edit.bed -o 5x_20k_90 --coverage 5 --identity_min 90 --length_mean 20000 --read_type pacbio --threads 40
#### 10X Coverage
VISOR LASeR -g ../resource/human_hs37d5.fasta -s chm1/ -b ../resource/SV_evaluation/SHORtS.LASeR.edit.bed -o 10x_20k_90 --coverage 10 --identity_min 90 --length_mean 20000 --read_type pacbio --threads 40
# 20X Coverage
VISOR LASeR -g ../resource/human_hs37d5.fasta -s chm1/ -b ../resource/SV_evaluation/SHORtS.LASeR.edit.bed -o 20x_20k_90 --coverage 20 --identity_min 90 --length_mean 20000 --read_type pacbio --threads 40
```

#### Modified. The commas (',') in all the read names were removed (replaced with vertical lines ('|')) for compatibility with NanoVar.
```
cd 3x_20k_90
samtools view -h sim.srt.bam | perl -lne '@row = split /\t/; $row[0] =~ s/,/|/g; print join ( "\t", @row );' | samtools view -Sb - -o sim.srt.edit.bam
cd 5x_20k_90
samtools view -h sim.srt.bam | perl -lne '@row = split /\t/; $row[0] =~ s/,/|/g; print join ( "\t", @row );' | samtools view -Sb - -o sim.srt.edit.bam
cd 10x_20k_90
samtools view -h sim.srt.bam | perl -lne '@row = split /\t/; $row[0] =~ s/,/|/g; print join ( "\t", @row );' | samtools view -Sb - -o sim.srt.edit.bam
cd 20x_20k_90
samtools view -h sim.srt.bam | perl -lne '@row = split /\t/; $row[0] =~ s/,/|/g; print join ( "\t", @row );' | samtools view -Sb - -o sim.srt.edit.bam
```

### Run NanoVar
```
cd ../3x_20k_90
nanovar sim.srt.edit.bam ../../resource/human_hs37d5.fasta nanovar -x pacbio-clr -t 24 --debug

cd ../5x_20k_90
nanovar sim.srt.edit.bam ../../resource/human_hs37d5.fasta nanovar -x pacbio-clr -t 24 --debug

cd ../10x_20k_90
nanovar sim.srt.edit.bam ../../resource/human_hs37d5.fasta nanovar -x pacbio-clr -t 24 --debug

cd ../20x_20k_90
nanovar sim.srt.edit.bam ../../resource/human_hs37d5.fasta nanovar -x pacbio-clr -t 24 --debug
```

### Run Sniffles
```
cd ../3x_20k_90
mkdir sniffles
cd sniffles
sniffles -m ../sim.srt.bam -v sniffles.vcf -s 1 -l 30 --genotype

cd ../../5x_20k_90
mkdir sniffles
cd sniffles
sniffles -m ../sim.srt.bam -v sniffles.vcf -s 2 -l 30 --genotype

cd ../../10x_20k_90
mkdir sniffles
cd sniffles
sniffles -m ../sim.srt.bam -v sniffles.vcf -s 3 -l 30 --genotype

cd ../../20x_20k_90
mkdir sniffles
cd sniffles
sniffles -m ../sim.srt.bam -v sniffles.vcf -s 4 -l 30 --genotype
```


### Process NanoVar VCFs and run Truvari (Modified)
```
cd ../../3x_20k_90/nanovar
# Modified. The '>' symbol was removed from the SVLEN field and entries with 'SVLEN=.' were removed due to incompatibility with Truvari.
perl -pe 's/SVLEN=>/SVLEN=/g' sim.srt.edit.nanovar.pass.vcf | grep 'SVLEN=\.' -v > sim.srt.edit.nanovar.pass.edit.vcf
bgzip -c sim.srt.edit.nanovar.pass.edit.vcf > sim.srt.edit.nanovar.pass.edit.vcf.gz
tabix sim.srt.edit.nanovar.pass.edit.vcf.gz
truvari bench -b ../../../resource/SV_evaluation/TOTAL.chm1.vcf.gz -c sim.srt.edit.nanovar.pass.edit.vcf.gz -o all -p 0 --sizemax 100000000

cd ../../5x_20k_90/nanovar
# Modified. The '>' symbol was removed from the SVLEN field and entries with 'SVLEN=.' were removed due to incompatibility with Truvari.
perl -pe 's/SVLEN=>/SVLEN=/g' sim.srt.edit.nanovar.pass.vcf | grep 'SVLEN=\.' -v > sim.srt.edit.nanovar.pass.edit.vcf
bgzip -c sim.srt.edit.nanovar.pass.edit.vcf > sim.srt.edit.nanovar.pass.edit.vcf.gz
tabix sim.srt.edit.nanovar.pass.edit.vcf.gz
truvari bench -b ../../../resource/SV_evaluation/TOTAL.chm1.vcf.gz -c sim.srt.edit.nanovar.pass.edit.vcf.gz -o all -p 0 --sizemax 100000000

cd ../../10x_20k_90/nanovar
# Modified. The '>' symbol was removed from the SVLEN field and entries with 'SVLEN=.' were removed due to incompatibility with Truvari.
perl -pe 's/SVLEN=>/SVLEN=/g' sim.srt.edit.nanovar.pass.vcf | grep 'SVLEN=\.' -v > sim.srt.edit.nanovar.pass.edit.vcf
bgzip -c sim.srt.edit.nanovar.pass.edit.vcf > sim.srt.edit.nanovar.pass.edit.vcf.gz
tabix sim.srt.edit.nanovar.pass.edit.vcf.gz
truvari bench -b ../../../resource/SV_evaluation/TOTAL.chm1.vcf.gz -c sim.srt.edit.nanovar.pass.edit.vcf.gz -o all -p 0 --sizemax 100000000

cd ../../20x_20k_90/nanovar
# Modified. The '>' symbol was removed from the SVLEN field and entries with 'SVLEN=.' were removed due to incompatibility with Truvari.
perl -pe 's/SVLEN=>/SVLEN=/g' sim.srt.edit.nanovar.pass.vcf | grep 'SVLEN=\.' -v > sim.srt.edit.nanovar.pass.edit.vcf
bgzip -c sim.srt.edit.nanovar.pass.edit.vcf > sim.srt.edit.nanovar.pass.edit.vcf.gz
tabix sim.srt.edit.nanovar.pass.edit.vcf.gz
truvari bench -b ../../../resource/SV_evaluation/TOTAL.chm1.vcf.gz -c sim.srt.edit.nanovar.pass.edit.vcf.gz -o all -p 0 --sizemax 100000000
```

### Process Sniffles VCFs and run Truvari (Modified)
```
cd ../../3x_20k_90/sniffles
grep '#' sniffles.vcf > head
grep -v '#' sniffles.vcf > body
sort -k 1,1 -k 2,2n body > Body
cat head Body > sniffles.vcf
rm head Body body
bgzip -c sniffles.vcf > sniffles.vcf.gz
tabix sniffles.vcf.gz
truvari bench -b ../../../resource/SV_evaluation/TOTAL.chm1.vcf.gz -c sniffles.vcf.gz -o all -p 0 --sizemax 100000000

cd ../../5x_20k_90/sniffles
grep '#' sniffles.vcf > head
grep -v '#' sniffles.vcf > body
sort -k 1,1 -k 2,2n body > Body
cat head Body > sniffles.vcf
rm head Body body
bgzip -c sniffles.vcf > sniffles.vcf.gz
tabix sniffles.vcf.gz
truvari bench -b ../../../resource/SV_evaluation/TOTAL.chm1.vcf.gz -c sniffles.vcf.gz -o all -p 0 --sizemax 100000000

cd ../../10x_20k_90/sniffles
mv sniffles.vcf sniffles.original.vcf
# Modified. Entries with 'STRANDBIAS' in the 'FILTER' column were removed due to incompatibility with Truvari.
grep 'STRANDBIAS' -v sniffles.original.vcf > sniffles.vcf
grep '#' sniffles.vcf > head
grep -v '#' sniffles.vcf > body
sort -k 1,1 -k 2,2n body > Body
cat head Body > sniffles.vcf
rm head Body body
bgzip -c sniffles.vcf > sniffles.vcf.gz
tabix sniffles.vcf.gz
truvari bench -b ../../../resource/SV_evaluation/TOTAL.chm1.vcf.gz -c sniffles.vcf.gz -o all -p 0 --sizemax 100000000

cd ../../20x_20k_90/sniffles
mv sniffles.vcf sniffles.original.vcf
# Modified. Entries with 'STRANDBIAS' in the 'FILTER' column were removed due to incompatibility with Truvari.
grep 'STRANDBIAS' -v sniffles.original.vcf > sniffles.vcf
grep '#' sniffles.vcf > head
grep -v '#' sniffles.vcf > body
sort -k 1,1 -k 2,2n body > Body
cat head Body > sniffles.vcf
rm head Body body
bgzip -c sniffles.vcf > sniffles.vcf.gz
tabix sniffles.vcf.gz
truvari bench -b ../../../resource/SV_evaluation/TOTAL.chm1.vcf.gz -c sniffles.vcf.gz -o all -p 0 --sizemax 100000000
```

### The summary files of all Truvari analysis can be found [here](https://github.com/cytham/nv_benchmark_jiang/tree/main/truvari_summary).
