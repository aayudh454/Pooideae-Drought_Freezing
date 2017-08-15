# RNA seq of B. distachyon and N. pulchra with drought and freezing stress

### Author: Aayudh Das

Login info: **ssh aadas@bluemoon-user2.uvm.edu**

### Data structure (*B. distachyon and N. pulchra* data)

**Treatment:** Both B. distachyon and N. pulchra has three replicates for "control" (BC1-BC3 and NC1-NC3), three for drought (BD1-BD3 and ND1-ND3) three for freezing (BF1-BF4, and NF1-NF3). Wonderful data set.

**Method for sequencing:** Single eneded and stranded sequencing.

## Table of contents    
* [Page 1: 2017-08-05](#id-section1). Moving files and trimming by Trimmomatic-0.36

* [Page 2: 2017-08-05](#id-section2). Trimming for all files and Concatenation

* [Page 3 2017-08-06](#id-section3). Assembly by Trinity 2.1.1

* [Page 4 2017-08-07](#id-section4). Transcript quantification by RSEM 

* [Page 5 2017-08-08](#id-section5). Build Transcript and Gene Expression Matrices

* [Page 6 2017-08-08](#id-section6). Differential Expression Analysis (DeSeq2)

* [Page 7 2017-04-14](#id-section7). Gene Annotation with blastP

* [Page 8 2017-04-25](#id-section8). Principal Component Analysis (PCA)

* [Page 9 2017-05-02](#id-section9). 

  ​

------
<div id='id-section1'/>
###Page 1: 2017-08-05. Moving files and trimming by Trimmomatic-0.36

#### 1. Where are our files?

After login you should **cd** space the location path. **ll** to see what are the files present

```
[aadas@bluemoon-user2 ~]$ cd /gpfs2/scratch/djshirle/MPS/170728_SNL128_0160_AHMTLCBCXY/samples_out
[aadas@bluemoon-user2 samples_out]$ ll
-rw-r--r-- 1 aadas usr 6106 Aug  1 16:37 170729_demux_sheet-Preston.csv.quality_table_by_lane.txt
-rw-r--r-- 1 aadas usr 3368 Aug  1 16:37 170729_demux_sheet-Preston.csv.quality_table_by_sample.txt
drwxr-xr-x 3 aadas usr  512 Aug  1 16:37 BC1
drwxr-xr-x 3 aadas usr  512 Aug  1 16:37 BC2
drwxr-xr-x 3 aadas usr  512 Aug  1 16:37 BC3
drwxr-xr-x 3 aadas usr  512 Aug  1 16:37 BD1
drwxr-xr-x 3 aadas usr  512 Aug  1 16:37 BD2
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 BD3
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 BF1
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 BF2
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 BF4
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 NC1
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 NC2
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 NC3
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 ND1
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 ND2
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 ND3
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 NF1
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 NF2
drwxr-xr-x 3 aadas usr  512 Aug  1 16:38 NF3
```

#### 2. How can you move all the files to your own directory from the server?

```
[aadas@bluemoon-user2 npulbrachyaug2017]$ cp -r /gpfs2/scratch/djshirle/MPS/170728_SNL128_0160_AHMTLCBCXY/samples_out/* . &
```

cp=copy, *= means everything (for individual file just write the name of the file instead of *)

#### 3. Let's see what files are present in BC1 folder "Control" (Brachy) 

```
[aadas@bluemoon-user2 BC1]$ ll
total 2200432
drwxr-xr-x 4 aadas usr        512 Aug  1 16:37 jcpresto_PooidFreezeDrought_20170614_BC1_R1_fastqc
-rw-r--r-- 1 aadas usr     218646 Aug  1 16:37 jcpresto_PooidFreezeDrought_20170614_BC1_R1_fastqc.zip
-rw-r--r-- 1 aadas usr 2253018076 Aug  1 16:37 jcpresto_PooidFreezeDrought_20170614_BC1_R1.fastq.gz
-rw-r--r-- 1 aadas usr        165 Aug  1 16:37 jcpresto_PooidFreezeDrought_20170614_BC1_R1.fastq.gz.md5sum

[aadas@bluemoon-user2 ~]$ cd nassellaBrachy_drought_freezing/trimming_nassellaBrachy/
```

#### 4. Copy fastq.gz files of all the replicates into trimming folder 

Now you need **R1.fastq.gz** file, so copy them to your **trimming_nassellaBrachy** folder. So, finish copying all the files.

```
[aadas@bluemoon-user2 npulbrachyaug2017]$ cd BC2/
[aadas@bluemoon-user2 BC2]$ cp jcpresto_PooidFreezeDrought_20170614_BC2_R1.fastq.gz ~/nassellaBrachy_drought_freezing/trimming_nassellaBrachy/
```

#### 5. Now create a script for running the trimming program

you should be in your **trimming** folder. Now to create a blank script

```
[aadas@bluemoon-user2 trimming_nassellaBrachy]$ vi trim.sh
```

Now press **i** to insert 

Copy and paste everything present in the new script from this sample one

and **esc**, **:wq** to save and quit. 

**:q!** only quit without saving.

look how I edited the things

```
#!/bin/bash

######## This job needs 1 nodes, 2 processors total
#PBS -q poolmemq
#PBS -l nodes=1:ppn=2,mem=16gb,vmem=18gb
# it needs to run for 6 hours
#PBS -l walltime=06:00:00
#PBS -N job
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea
###LOAD JAVA MODULE AVAILABLE FROM THE CLUSTER, YOU MAY WANT TO CHECK FIRST
module load java-sdk/sun-jdk-1.6.0.12
ulimit -s unlimited
###CHANGE THE DIRECTORY ACCORDINGLY, THE FOLLOWING SETTINGS ARE FOR MY ACCOUNT
SOFTWARE=/users/a/a/aadas/Bin/Trimmomatic-0.36
workDIR=/users/a/a/aadas/nassellaBrachy_drought_freezing/trimming_nassellaBrachy
cd $workDIR
#####TRIMMING COMMANDS AND PARAMETERS
java -jar $SOFTWARE/trimmomatic-0.36.jar SE -phred33 $workDIR/jcpresto_PooidFreezeDrought_20170614_NF3_R1.fastq.gz $workDIR/NF3_trimmed.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:40
```

#### What does all this things actually mean?

**Based on the sequencing** (paired end or signle end)  you need to follow this protocol BUT modify other factors! Don't copy deto!!!!!

#### Paired End:

`java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

This will perform the following:

- Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
- Remove leading low quality or N bases (below quality 3) (LEADING:3)
- Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
- Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
- Drop reads below the 36 bases long (MINLEN:36)

#### Single End:

`java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

This will perform the same steps, using the single-ended adapter file

#### 6. Now submit the job and other editing

i. You need to check **Submitting Jobs to the Cluster**-

https://www.uvm.edu/~vacc/?Page=userguide.php

first part of the script explains that details

ii. Now you need to specify where your software is present i.e. the Trimmomatic-0.36 which is in your main directory **/users/a/a/aadas/Bin/Trimmomatic-0.36**

iii. Specify your working directory-

Just **pwd** where ever your files are and then copy that in the script.

**/users/a/a/aadas/nassellaBrachy_drought_freezing/trimming_nassellaBrachy**

iv. TRIMMING COMMANDS AND PARAMETERS

a. change software version from as **trimmomatic-0.36**

b. Now the **first file is your input file**, so after $workDIR/"name of the file" space

c. Last file is your **output file**.  

**$workDIR/NF3_trimmed.fq.gz**

#### 7. Make your script executable

you should be the folder where you saved your script

```
[aadas@bluemoon-user2 Ba]$ chmod 700 trimmomatic.sh
```

**700**=file's owner may read, write, and execute the file.

#### 8. Submit your job and check status of your job

```
[aadas@bluemoon-user2 Ba]$ qsub trimmomatic.sh 
```

check status

```
[aadas@bluemoon-user2 Ba]$ showq -u aadas
```

all jobs

```
[aadas@bluemoon-user2 Ba]$ showq
```
hang on there it's gonna take 6hrs. You can submit multiple jobs.

-----
<div id='id-section2'/>
### Page 2: 2017-08-05. Trimming for all files and Concatenation 

###### Some moving tips 

**Rename a file**

```
mv oldname newname
mv *.trimmo.fq.gz ../Brachyleytrum_aristosum
```

#### 1. Moving all files to 'trimming_nassellaBrachy' folder, edit the script and execute that

```
[aadas@bluemoon-user2 trimming_nassellaBrachy]$ qsub trim.sh
[aadas@bluemoon-user2 trimming_nassellaBrachy]$ showq -u aadas
```

#### Copy only the trimmo.fq.gz files to a new folder for assembly

```
[aadas@bluemoon-user2 trimming_nassellaBrachy]$ cp *.trimmo.fq.gz ../Brachyleytrum_aristosum
```

#### 2. Execute concatenation 

```
[aadas@bluemoon-user2 trimming_nassellaBrachy]$ zcat *trimmed.fq.gz > brachyNpul.R1.trimmed.fq &
```

It takes about 30 min. Don't think it's done right away!

#### 3. After finishing the concatenation check the number of reads

```
[aadas@bluemoon-user2 trimming_nassellaBrachy]$ grep -c "@" brachyNpul.R1.trimmed.fq 
```

It shows **322498031** reads. 

Like I compared my previous analysis of Brachyletrium, it turns out that 72805944/168097158 is only **43%** .

If it's too low then modify LEADING:20 TRAILING:20 accordingly to LEADING:3 TRAILING:3

check that how much you get if you play with the parameters.

#### 4. Move the catenated file to a new folder where you will do assembly

``` 
[aadas@bluemoon-user2 trimming_nassellaBrachy]$ cp brachyNpul_catenated.fq ~/nassellaBrachy_drought_freezing/assembly_npulBdis/
```

------

<div id='id-section3'/>

### Page 3: 2017-08-06. Assembly by Trinity 2.1.1

### Some tips

to delete a line: **esc** and **d+d**

https://www.cs.colostate.edu/helpdocs/vi.html

Save all the commands as a history.txt

```
[aadas@bluemoon-user2 ~]$ history > history.txt
```

##### When you download something in mac and upload that to server

```
scp bowtie-1.2-linux-x86_64.zip aadas@bluemoon-user2.uvm.edu:~/
```

##### Delete an entire directory

```
rm -rf filename
```

#### 1. Install trinity (using Trinity 2.1.1) 

Go to - https://github.com/trinityrnaseq/trinityrnaseq/releases and then copy link address of **Source code (tar.gz)**

```
[aadas@bluemoon-user2 ~]$ wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.tar.gz
```

Now create a new directory and move Trinity to that

```
[aadas@bluemoon-user2 ~]$ mkdir Bin
[aadas@bluemoon-user2 ~]$ mv Trinity-v2.4.0.tar.gz Bin/
```

unzip the tar file

```
[aadas@bluemoon-user2 Bin]$ tar -zxvf Trinity-v2.4.0.tar.gz
```

- -z: Compress the archive with g**z**ip.
- x ; decompress or extraction
- -v: Display progress in the terminal while creating the archive, also known as “**v**erbose” mode. The v is always optional in these commands, but it’s helpful.
- -f: Allows you to specify the **f**ilename of the archive.

now install 

```
[aadas@bluemoon-user2 Bin]$ cd trinityrnaseq-Trinity-v2.4.0/
[aadas@bluemoon-user2 trinityrnaseq-Trinity-v2.4.0]$ make 
```

but when downloading Trinity v2.0.6

````
[aadas@bluemoon-user2 trinityrnaseq-2.0.6]$ make plugins
````

edit the script (change WORKINGDIR file names). Here you need to include the module load but first you need to check what modules are availble in the cluster.

```
[aadas@bluemoon-user2 Brachyleytrum_aristosum]$ module avail
```

Then include the samtools and bowtie version, so that your script can use that program.

Increase the cpu load to 16 in **ppn=16** and also in the end **CPU 16** 

#### 2. Sript for assembly for single ended by trinity 2.1.1

```
#!/bin/bash

#PBS -l nodes=1:ppn=24,mem=256G,vmem=288G
#PBS -q poolmemq
# it needs to run for 6 hours
#PBS -l walltime=30:00:00
#PBS -N trinityv211
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea
module load samtools-1.3.1-gcc-6.3.0-e5jw5u4
module load bowtie2-2.2.5-gcc-6.3.0-daskah5
export PATH="/users/a/a/aadas/Bin/bowtie-1.1.1:$PATH"
#ulimit -s unlimited

SOFTWAREDIR=/users/a/a/aadas/Bin/trinityrnaseq-2.1.1
WORKINGDIR=/users/a/a/aadas/nassellaBrachy_drought_freezing/assembly_npulBdis
cd $WORKINGDIR

/users/a/a/aadas/Bin/trinityrnaseq-2.1.1/Trinity --seqType fq --normalize_reads --max_memory 256G --single /users/a/a/aadas/nassellaBrachy_drought_freezing/assembly_npulBdis/brachyNpul_catenated.fq --CPU 24
```

##### What all those mean?

- Use all reads from an individual (all conditions) to capture most genes
- Read files may be gzipped (as in this example) or not (then they should not have the “.gz” ending)
- Paired-end reads specified with **--left** and -**-right**. If only single-end, use **--single** instead.
- 256G is the maximum memory to be used at any stage which allows memory limitation (jellyfish, sorting, etc.)
- At most 24 CPU cores will be used in any stage.

#### First time before you execute the job with a new Trinity version

```
  936  2017-04-05 10:28:20 module load samtools-1.3.1-gcc-6.3.0-e5jw5u4
  937  2017-04-05 10:28:32 export PATH="/users/a/a/aadas/Bin/bowtie-1.2:$PATH"
  938  2017-04-05 10:28:40 export PATH="/users/a/a/aadas/Bin/jre1.7.0_80/bin:$PATH"
  939  2017-04-05 10:28:46 export PATH="/users/a/a/aadas/Bin/jre1.7.0_80/bin/java:$PATH"
  940  2017-04-05 10:29:10 cd ..
  941  2017-04-05 10:29:19 cd Bin/
  942  2017-04-05 10:29:20 ll
  943  2017-04-05 10:29:32 cd trinityrnaseq-Trinity-v2.3.2
  944  2017-04-05 10:29:35 make
  945  2017-04-05 10:34:06 cd ~/
```

Submit the job and view. If you see that 30% completed after 30hrs then qsub again and it will catch up from where it was left before.

#### 3. Check the assembly file

```
# number of sequence present
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta | less
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | less
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | sort -u | less
[aadas@bluemoon-user2 assembly_npulBdis]$ grep ">" brachyNpul_Trinity211.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | sort -u | wc -l
```

**check no. of seq**

```
[aadas@bluemoon-user2 assembly_npulBdis]$ grep ">" npulBdis_Trinity211.fasta -c
```

**423850**

Check back after your job is done

```
[aadas@bluemoon-user2 assembly_npulBdis]$ grep ">" npulBdis_Trinity211.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | sort -u | wc -l
```

**274312**

------

<div id='id-section4'/>

### Page 4: 2017-08-07. Transcript quantification by RSEM

RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome

Script for transcript quantification.

* Make sure your **Trinity.fasta** file in the directory
* Make sure that **.trimmed.fq.gz** is also in the FASTQ_DIR

### For Brachypodium_control

```
#!/bin/bash

######## This job needs 1 nodes, 4 processors total
#PBS -l nodes=1:ppn=4,pmem=8gb,pvmem=9gb
#PBS -l walltime=30:00:00
#PBS -N brachypodium_control
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea
module load samtools-1.3.1-gcc-6.3.0-e5jw5u4

TRINITY_HOME=/users/a/a/aadas/Bin/trinityrnaseq-2.1.1

export PATH=/users/a/a/aadas/Bin/bowtie-1.1.1:$PATH
export PATH=/users/a/a/aadas/Bin/RSEM-1.2.19:$PATH

FASTQ_DIR=/users/a/a/aadas/nassellaBrachy_drought_freezing/rsem_npulBdis
cd $FASTQ_DIR

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts npulBdis_Trinity211.fasta --seqType fq --single $FASTQ_DIR/BC1_trimmed.fq.gz --est_method RSEM --aln_method bowtie --thread_count 4 --SS_lib_type R --trinity_mode --prep_reference --output_dir $FASTQ_DIR --output_prefix Brachypodium_control01

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts npulBdis_Trinity211.fasta --seqType fq --single $FASTQ_DIR/BC2_trimmed.fq.gz --est_method RSEM --aln_method bowtie --thread_count 4 --SS_lib_type R --trinity_mode --prep_reference --output_dir $FASTQ_DIR --output_prefix Brachypodium_control02

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts npulBdis_Trinity211.fasta --seqType fq --single $FASTQ_DIR/BC3_trimmed.fq.gz --est_method RSEM --aln_method bowtie --thread_count 4 --SS_lib_type R --trinity_mode --prep_reference --output_dir $FASTQ_DIR --output_prefix Brachypodium_control03
```

submit the job

``` 
[aadas@bluemoon-user2 rsem_npulBdis]$ qsub brachypodium_control.sh 
```

DO SAME script for all other treatments and also for Nassella.

#### Open R in UNIX

You should be in the directory where your files seats.

```
[aadas@bluemoon-user2 rsem_npulBdis]$ module load r-3.3.2-gcc-6.3.0-bmdvb4s
[aadas@bluemoon-user2 rsem_npulBdis]$ R
```

few tips: Type 'q()' to quit R.

```
> q()
Save workspace image? [y/n/c]: n
```

The above `R` commands list the top **10 highest expressed genes** (shown below).

```
> data = read.table("Brachypodium_control01.genes.results", header=T, stringsAsFactors=F)
> idx = order(data[,"TPM"], decreasing=T)
> data[idx[1:10], c("gene_id", "expected_count", "TPM")]
                     gene_id expected_count       TPM
208017 TRINITY_DN75056_c7_g1     7959401.35 215296.86
212426 TRINITY_DN75795_c6_g2     8291434.79 185797.31
226295 TRINITY_DN77927_c4_g3       72412.00   9846.73
249661 TRINITY_DN81029_c8_g1       50093.38   8856.78
23512   TRINITY_DN1199_c0_g1        6215.77   8454.08
108931 TRINITY_DN35614_c0_g1        8702.09   5739.32
224027 TRINITY_DN77596_c6_g6       23774.50   5453.07
188162 TRINITY_DN71287_c3_g1       31410.34   5402.16
251927 TRINITY_DN81269_c3_g1      211601.00   5217.94
120536 TRINITY_DN42492_c0_g1        3860.62   4931.95
```

Do this for all replicates!

------

<div id='id-section5'/>

### Page 5: 2017-08-08. Build Transcript and Gene Expression Matrices

Terms:

**<u>FPKM</u>: fragments per kilobase transcript length per million fragments mapped**

**<u>TPM</u>: transcripts per million transcripts**

Using the transcript and gene-level abundance estimates for each of your samples, construct a matrix of counts and a matrix of normalized expression values using the following script:

#### For genes of *Brachypodium* control, drought and freezing

```
[aadas@bluemoon-user2 rsem_npulBdis]$ ~/Bin/trinityrnaseq-2.1.1/util/abundance_estimates_to_matrix.pl --est_method RSEM Brachypodium_control01.genes.results Brachypodium_control02.genes.results Brachypodium_control03.genes.results Brachypodium_drought01.genes.results Brachypodium_drought02.genes.results Brachypodium_drought03.genes.results Brachypodium_freezing01.genes.results Brachypodium_freezing02.genes.results Brachypodium_freezing03.genes.results --out_prefix Brachypodium.genes 
```

#### For genes of *Nassella* control, drought and freezing

```
~/Bin/trinityrnaseq-2.1.1/util/abundance_estimates_to_matrix.pl --est_method RSEM Nassella_control01.genes.results Nassella_control02.genes.results Nassella_control03.genes.results Nassella_droughtl01.genes.results Nassella_drought02.genes.results Nassella_drought03.genes.results Nassella_freezing01.genes.results Nassella_freezing02.genes.results Nassella_freezing03.genes.results --out_prefix Nassella.genes 
```

#### For genes of *Brachypodium* and *Nassella* control, drought and freezing

```
[aadas@bluemoon-user2 rsem_npulBdis]$ ~/Bin/trinityrnaseq-2.1.1/util/abundance_estimates_to_matrix.pl --est_method RSEM Brachypodium_control01.genes.results Brachypodium_control02.genes.results Brachypodium_control03.genes.results Brachypodium_drought01.genes.results Brachypodium_drought02.genes.results Brachypodium_drought03.genes.results Brachypodium_freezing01.genes.results Brachypodium_freezing02.genes.results Brachypodium_freezing03.genes.results Nassella_control01.genes.results Nassella_control02.genes.results Nassella_control03.genes.results Nassella_droughtl01.genes.results Nassella_drought02.genes.results Nassella_drought03.genes.results Nassella_freezing01.genes.results Nassella_freezing02.genes.results Nassella_freezing03.genes.results --out_prefix Brachypodium_Nassella.genes 
```

#### Now smilar way you can do this for isoforms

```
~/Bin/trinityrnaseq-2.1.1/util/abundance_estimates_to_matrix.pl --est_method RSEM Brachyleytrum_precold01.isoforms.results Brachyleytrum_precold02.isoforms.results Brachyleytrum_precold03.isoforms.results Brachyleytrum_coldshock01.isoforms.results Brachyleytrum_coldshock02.isoforms.results Brachyleytrum_coldshock03.isoforms.results Brachyleytrum_sixweekcold01.isoforms.results Brachyleytrum_sixweekcold02.isoforms.results --out_prefix Brachyleytrum.isoforms
```

#### Counting Numbers of Expressed Transcripts or Genes

```
~/Bin/trinityrnaseq-2.1.1/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \ Brachypodium_Nassella.genes.TPM.not_cross_norm | tee Brachypodium_Nassella.genes.TPM.not_cross_norm.counts_by_min_TPM

~/Bin/trinityrnaseq-2.1.1/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
          trans_matrix.TPM.not_cross_norm | tee trans_matrix.TPM.not_cross_norm.counts_by_min_TPM
```

The above table indicates that we have 52,624 'genes' that are expressed by at least 1 TPM in any one of the many samples in this expression matrix.

Plotting the number of 'genes' (or 'transcripts') as a function of minimum TPM threshold, we can see that the vast majority of all expressed features have very little expression support. Using R (or your own favorite data analysis package), we might extrapolate the number of expressed 'genes' based on the trend prior to the massive influx of lowly expressed transcripts:

```
setwd("~/Desktop")
list.files()
data = read.table("Brachypodium_Nassella.genes.TPM.not_cross_norm.counts_by_min_TPM", header=T)
plot(data, xlim=c(-100,0), ylim=c(0,100000), t='b')
```

```
# extract the data between 10 TPM and 100 TPM
filt_data = data[data[,1] > -100 & data[,1] < -10,] 
# perform a linear regression on this filtered subset of the data
fit = lm(filt_data[,2] ~ filt_data[,1])
print(fit)

Call:
  lm(formula = filt_data[, 2] ~ filt_data[, 1])

Coefficients:
  (Intercept)  filt_data[, 1]  
19361.4           195.5 

# add the linear regression line to the plot 
abline(fit, col='green', lwd=3)
```

The linear regression allows us to extrapolate (based on the Y-intercept) that we have 19361 'genes', which is a far better guess than our count of 52,624 'genes' having at least 1 TPM in any sample, and certainly better than the 1.4 million 'genes' that were assembled. 

------

<div id='id-section6'/>

### Page 6: 2017-08-08. Differential Expression Analysis (DeSeq2)

First lets get 'R' working

```
[aadas@bluemoon-user2 ~]$ module load r-3.3.2-gcc-6.3.0-bmdvb4s
```

```
R
 > source("http://bioconductor.org/biocLite.R")
 > biocLite('edgeR')
 > biocLite('limma')
 > biocLite('DESeq2')
 > biocLite('ctc')
 > biocLite('Biobase')
 > install.packages('gplots')
 > install.packages('ape')
```

### 1. Identifying DE features: With biological replicates (PREFERRED)

Be sure to have a tab-delimited 'samples_described.txt' file that describes the relationship between samples and replicates. For example:

just create a text file and copy paste this

```
conditionA   Brachypodium_control01
conditionA   Brachypodium_control02
conditionA   Brachypodium_control03

conditionB   Brachypodium_drought01
conditionB   Brachypodium_drought02
conditionB   Brachypodium_drought03

conditionC   Brachypodium_freezing01
conditionC   Brachypodium_freezing02
conditionC   Brachypodium_freezing03

conditionD   Nassella_control01
conditionD   Nassella_control02
conditionD   Nassella_control03

conditionE   Nassella_droughtl01
conditionE   Nassella_drought02
conditionE   Nassella_drought03

conditionF   Nassella_freezing01
conditionF   Nassella_freezing02
conditionF   Nassella_freezing03
```

#### Voom method for differential gebe expression

Any of the available methods support analyses containing biological replicates. Here, for example, we again choose voom within the limma package.

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Brachypodium_Nassella.genes.counts.matrix --method voom --samples_file samples_described.txt  
```

### 2. Extracting and clustering differentially expressed transcripts (HEATMAP)

An initial step in analyzing differential expression is to extract those transcripts that are most differentially expressed (most significant FDR and fold-changes) and to cluster the transcripts according to their patterns of differential expression across the samples. You should do this **inside the voom folder**. 

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix Brachypodium_Nassella.genes.counts.matrix -P 1e-3 -C 2 --samples samples_described.txt
```

you might need to update the q value in library.

which will extract all genes that have P-values at most 1e-3 and are at least 2^2 fold differentially expressed. For each of the earlier pairwise DE comparisons, this step will generate the following files:

if it's more than 10000 diff genes-

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix Brachypodium_Nassella.genes.counts.matrix -P 1e-3 -C 2 --samples samples_described.txt --max_genes_clust 50000
```

### 3. Automatically Partitioning Genes into Expression Clusters

There are three different methods for partitioning genes into clusters:

- use K-means clustering to define K gene sets. (use the -K parameter). This does not leverage the already hierarchically clustered genes as shown in the heatmap, and instead uses a least-sum-of-squares method to define exactly k gene clusters.
- cut the hierarchically clustered genes (as shown in the heatmap) into exactly K clusters.
- (Recommended) cut the hierarchically clustered gene tree at --Ptree percent height of the tree.

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl \
                                    -R diffExpr.P1e-6_C7.1.matrix.RData --Ptree 60
```

You should do this **inside the voom folder**. 

### 4. Principal Component Analysis (PCA)

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/PtR --matrix Brachypodium_Nassella.genes.counts.matrix \
    -s samples_described.txt --log2 --prin_comp 3
```

The --prin_comp 3 indicates that the first three principal components will be plotted, as shown above, with PC1 vs. PC2 and PC2 vs. PC3. In this example, the replicates cluster tightly according to sample type, which is very reassuring.

If you have replicates that are clear outliers, you might consider removing them from your study as potential confounders. If it's clear that you have a [batch effect](http://www.nature.com/nrg/journal/v11/n10/full/nrg2825.html), you'll want to eliminate the batch effect during your downstream analysis of differential expression.

### Now only Brachypodium

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Brachypodium.genes.counts.matrix --method voom --samples_file samples_described.txt  
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix Brachypodium.genes.counts.matrix -P 1e-2 -C 2 --samples samples_described.txt
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -R diffExpr.P1e-2_C2.matrix.RData --Ptree 60
```

### Now only Nassella

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Nassella.genes.counts.matrix --method voom --samples_file samples_described.txt  
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix Nassella.genes.counts.matrix -P 1e-3 -C 2 --samples samples_described.txt
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix Nassella.genes.counts.matrix -P 1e-2 -C 2 --samples samples_described.txt
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -R diffExpr.P1e-3_C2.matrix.RData --Ptree 60
```

#### Nassella vs Brachypodium control-drought

```
[aadas@bluemoon-user2 rsem_npulBdis]$ ~/Bin/trinityrnaseq-2.1.1/util/abundance_estimates_to_matrix.pl --est_method RSEM Brachypodium_control01.genes.results Brachypodium_control02.genes.results Brachypodium_control03.genes.results Brachypodium_drought01.genes.results Brachypodium_drought02.genes.results Brachypodium_drought03.genes.results  Nassella_control01.genes.results Nassella_control02.genes.results Nassella_control03.genes.results Nassella_droughtl01.genes.results Nassella_drought02.genes.results Nassella_drought03.genes.results --out_prefix BrachypodiumNassella_DroughtControl.genes 
```

just create a text file and copy paste this

```
conditionA   Brachypodium_control01
conditionA   Brachypodium_control02
conditionA   Brachypodium_control03

conditionB   Brachypodium_drought01
conditionB   Brachypodium_drought02
conditionB   Brachypodium_drought03

conditionD   Nassella_control01
conditionD   Nassella_control02
conditionD   Nassella_control03

conditionE   Nassella_droughtl01
conditionE   Nassella_drought02
conditionE   Nassella_drought03
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix BrachypodiumNassella_DroughtControl.genes.counts.matrix --method voom --samples_file samples_described.txt
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix BrachypodiumNassella_DroughtControl.genes.counts.matrix -P 1e-2 -C 2 --samples samples_described.txt
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix BrachypodiumNassella_DroughtControl.genes.counts.matrix -P 1e-3 -C 2 --samples samples_described.txt
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -R diffExpr.P1e-3_C2.matrix.RData --Ptree 60
```

### Extended version of Deseq analysis in R

```
setwd("~/Dropbox/Aayudh PhD/Brachy_Npul_transcriptomics")
list.files()
library("DESeq2")
library("ggplot2")

countsTable <- read.delim('Brachypodium_Nassella.genes.counts.matrix.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData=as.matrix(countsTable)
storage.mode(countData) = "integer"
head(countData)

conds <- read.delim("colsData_all.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)
dim(countData)
dim(colData)
#subset
colData_BrachyCvsD = colData[1:6,]
colData_BrachyCvsD
dim(countData)
dim(colData_BrachyCvsD)
countData_BrachyCvsD = countData[,1:6]
storage.mode(countData_BrachyCvsD) = "integer"
dim(countData_BrachyCvsD)
dim(colData_BrachyCvsD)
#model 
dds <- DESeqDataSetFromMatrix(countData = countData_BrachyCvsD,
                              colData = colData_BrachyCvsD,
                              design = ~ cond)
dds
dim(dds)
# Filtering to remove rows with 0 reads
dds <- dds[ rowSums(counts(dds)) > 1, ]
dim(dds)
dds <- DESeq(dds) 
res <- results(dds)
#sorts according to pvalue
res <- res[order(res$padj),]
#Subsetting based on pvalue
x.sub <- subset(res, padj < 0.01)
head(x.sub)
sub_data=as.data.frame(x.sub)
write.csv(sub_data, file = "Brachypodium_ControlvsDrought.csv", row.names = T, quote = F)

head(res)
summary(res)
#MA plot
plotMA(res, main="Brachypodium:Control vs Drought", ylim=c(-5,7), xlim=c(1e+02, 1e+04))
#dispersion
plotDispEsts( dds, ylim = c(1e-4, 1e+1), xlim = c(1e+1, 1e+5) )
hist( res$pvalue, breaks=20, col="grey" )
#rlog transform
rld <- rlog( dds )
head( assay(rld) )
par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )
sampleDists <- dist(t(assay(rld)))
sampleDists
#Heatmap
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment,
                                     rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
#gene clustering heatmap
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",margin=c(5,10),
           trace="none", dendrogram="column",
           offsetRow = 0.25, cexRow=0.6,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

heatmap.2(assay(rld)[ topVarGenes, ],Colv=FALSE,dendrogram="row",scale="row",
          col=cm.colors(256),trace="none", margin=c(5,10), 
          lwid=c(1.5,2.0),
          keysize=1, offsetRow = 0.25, cexRow=0.6)
#Plot counts
d <- plotCounts(dds, gene="TRINITY_DN139561_c0_g1", intgroup=(c("cond", "rep")), returnData=TRUE)
d
p <- ggplot(d, aes(x= cond, y=count, shape = cond, color=cond)) + 
  theme_minimal() + theme(text = element_text(size=20), 
                          panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p

#Effects of transformations on the variance(http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#plot-counts)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
rld <- rlog(dds, blind=FALSE)
meanSdPlot(assay(rld))
#PCA
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("cond"))
#PCA2
ramp <- 1:3/3
cols <- c(rgb(ramp, 0, 0),
          rgb(0, ramp, 0),
          rgb(0, 0, ramp),
          rgb(ramp, 0, ramp))
print( plotPCA( rld, intgroup = c( "cond"), col=cols ) )
## Normalizing using the method for an object of class"CountDataSet" 
dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)
## Performing estimation of dispersion parameter
dds.disp <- estimateDispersions(dds.norm)
plotDispEsts(dds.disp)


```



------

<div id='id-section8'/>

### Page 8: 2017-08-11. Gene Annotation with blastP

Follow this link-http://transdecoder.github.io/

The *TransDecoder* utility is run on a fasta file containing the target transcript sequences. The simplest usage is as follows: Install it first!

#### Install blast

go to-ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/  

Download **linux.tar.gz**

```
$ export PATH=$PATH:/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin
```

#### Install databases

- UniProt is another protein database that is highly currated and combines the Swiss-Prot, TrEMBL and PIR-PSD databases.  Each protein has a unique uniprot ID that is associated with a wealth of other information, including Gene Ontology GO categories. See the documentation for more information [http://www.uniprot.org/help/uniref](http://www.uniprot.org/help/uniref). The various databases can be downloaded from here using ftp or wget [http://www.ebi.ac.uk/uniprot/database/download.html](http://www.ebi.ac.uk/uniprot/database/download.html). 

```
[aadas@bluemoon-user2 uniprot]$ makeblastdb -in uniprot_sprot.pep -dbtype prot
```

Search the peptides for protein domains using Pfam. This requires [hmmer3](http://hmmer.janelia.org/) and [Pfam](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz) databases to be installed.

**Briefly, to compile HMMER from source:**

```
  % tar zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
  % cd hmmer-3.1b2
  % ./configure
  % make
  % make check
```

**Pfam-**

Uncompress and prepare the Pfam database for use with *hmmscan* like so:

```
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

```
export PATH=$PATH:/users/a/a/aadas/Bin/hmmer-3.1b2-linux-intel-x86_64/binaries

hmmscan --cpu 8 --domtblout pfam.domtblout /users/a/a/aadas/Bin/Pfam-A.hmm /users/a/a/aadas/annotation/long_ORF/longest_orfs.pep
```



### Step 1: extract the long open reading frames by Transdecoder

```
#!/bin/bash

#PBS -N out.transDecoder-step1
#PBS -l nodes=1:ppn=4,pmem=10G,pvmem=12g
#PBS -j oe
#PBS -l walltime=10:00:00
#PBS -M aadas@uvm.edu
#PBS -m bea

export PATH="/users/a/a/aadas/Bin/TransDecoder-3.0.1/transdecoder_plugins/cdhit:$PATH"
export PATH="/users/a/a/aadas/Bin/TransDecoder-3.0.1:$PATH"
export PATH="/users/a/a/aadas/Bin/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH"
export PATH="/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin:$PATH"

transDecoder_dir=/users/a/a/aadas/Bin/TransDecoder-3.0.1
INPUT_DIR=/users/a/a/aadas/nassellaBrachy_drought_freezing/blastP
cd $INPUT_DIR

$transDecoder_dir/TransDecoder.LongOrfs -t $INPUT_DIR/npulBdis_Trinity211.fasta
```



Step 2 

```
/users/a/a/aadas/Bin/TransDecoder-3.0.1/TransDecoder.Predict -t Brachyleytrum_trinityv211.fasta
```

#### Install blast

go to-ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/  

Download **linux.tar.gz**

```
$ export PATH=$PATH:/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin
```

**make the database**

```
[aadas@bluemoon-user2 uniprot]$ makeblastdb -in uniprot_sprot.pep -dbtype prot
```

**script**

This single line using the blastp command below will compare your transcript fasta file

(-query) to the already formatted uniref90 database (-db).

You can enter 'blastp --help' for a list of the parameters.

We choose the tab-delimited output format (6) and to only help the top hit (-max_target_seqs) and only if it has a minimum evalue of 0.001.

```
#!/bin/bash

blastp -query /users/a/a/aadas/annotation/long_ORF/longest_orfs.pep \
       -db /users/a/a/aadas/annotation/uniprot/uniprot_sprot.pep \
       -out ~/Brachy_vs_uniprot.outfmt6 \
       -outfmt 6 \
       -evalue 1e-3 \
       -max_target_seqs 1
```

now **bash** the script

### BlastP Search

If you want to submit as a job to VACC

```
#!/bin/bash

######## This job needs 1 nodes, 4 processors total
#PBS -l nodes=1:ppn=4,pmem=8gb,pvmem=9gb
#PBS -l walltime=30:00:00
#PBS -N blastp
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea

export PATH=$PATH:/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin

blastp -query /users/a/a/aadas/annotation/long_ORF/longest_orfs.pep \
       -db /users/a/a/aadas/annotation/uniprot/uniprot_sprot.pep \
       -out ~/Brachy_vs_uniprot.outfmt6 \
       -outfmt 6 \
       -evalue 1e-3 \
       -max_target_seqs 1
```

### Pfam Search

### Step 2 (run two scripts)

```
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


my $file;               #BLAST query sequences
my $size= 2000;         #Number of sequences per task
my $type;               #Blast program.
my $database;           #path to database
my $eval= 1e-5;         #BLAST e-value cutoff
my $outputFormat= 6;    #BLAST output format
my $outputDir="blast_out";      #output directory
my $help;
my $input_dir=`pwd`;
my $blast_dir="/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin";
my $hmm_dir="/users/a/a/aadas/Bin/hmmer-3.1b2-linux-intel-x86_64/binaries";
my $hmmDB_dir="/users/a/a/aadas/blastp_Brachyleytrum/database";
my $hmm="hmmscan";
my $hmmoutputDir="hmmscan_out";
GetOptions(
        'query=s'       => \$file,
        'num_seqs=i'    => \$size,
        'program=s'     => \$type,
        'database=s'    => \$database,
        'eval=f'        => \$eval,
        'm=i'           => \$outputFormat,
        'output_dir=s'  => \$outputDir,
#       'rcc_queue=s'   => \$queue,
#       'combine'       => \$autoCombine,
        'help'          => \$help
);

my $usage = <<__EOUSAGE__;

###################################################################################
#       Batch Blast: Task Array
###################################################################################
#
#  --query <string>             File containing query sequences
# 
#  --program <string>           The BLAST program to use
#
#  --database <string>          The location of the BLAST database
#
# Optional:
# 
#  --num_seqs <integer>         The number of sequences per data sub-set
#                               Default: 1000
#--eval <float>               The BLAST e-value cutoff
#                               Default: 1e-10
#
#  --m <integer>                BLAST output format (8 for tablular)
#                               Default: 8 - Tabular Output
#
#  --output_dir <string>        Output directory for BLAST results
#                               Default: blast_out
#
#  --combine                    Automatically combine and remove the split.*
#                               directories.
#                               Default: False
#
#  --job_name_prefix            A prefix to differentiate this batch_blast run
#                               from another.
#                               
#  --rcc_queue                  SGE queue to run the BLAST job on
#                               Default: rcc-30d
#
#  --help  
###################################################################################
#  
#  To pass other arguments to blastall use -- <args> AFTER all required
#  arguments.
#
###################################################################################

__EOUSAGE__
;

if(!$file || !$type || !$database || $help){
        die($usage);
}

my $seqid;
my $seq;

my $seq_counter=0;
my $split_count=1;

open my $infile, "<", $file;
while(<$infile>){
        chomp;
        #is this line a new sequence header?
        if(/^>/){
                #if there is a storred sequence then write it out to file first
                if($seqid){
                        mkdir "split.".$split_count;
                        chdir ("split.$split_count");
                        open my $OUT, ">>", "split.$split_count.fasta";
                        print $OUT "$seqid\n$seq\n";
                        close $OUT;
                        if($seq_counter == $size){
                                $seq_counter = 0;
                                $split_count++;

                        }
                        chdir ("../");
                }
                $seq=();
                $seqid = $_;
                $seq_counter++;
        }
        #if not then keep building the sequence
        else{
                $seq .= $_;
        }
}
}
#necessarily this loop exits before the last sequence is written. Write it now.
if($seqid){
        mkdir "split.".$split_count;
        chdir ("split.$split_count");
        open my $OUTsplit, ">>", "split.$split_count.fasta";
        print $OUTsplit "$seqid\n$seq\n";
        close $OUTsplit;
        if($seq_counter == $size){
                $seq_counter = 0;
                $split_count++;
        }
        chdir ("../");
}
close $infile;
#if the output directory doesn't exist then make it
if(! -e $outputDir){
        `mkdir $outputDir`;
}

#print the task-array script
for(my $i=1;$i<=$split_count;$i++){
     open my $SUB_SCRIPTS, ">", "$type-part-$i.sh" or die();
     print $SUB_SCRIPTS "#!/bin/bash\n",
                        "#PBS -N out.$type.part-$i\n",
                        "#PBS -l nodes=1:ppn=1,pmem=10G,pvmem=12g\n",
                        "#PBS -j oe\n",
                        "#PBS -l walltime=30:00:00\n",
                        "#PBS -M aadas\@uvm.edu\n",
                        "#PBS -m bea\n",
                        "\n",
                        "INPUT_DIR=$input_dir\n",
                        "BLAST_DIR=$blast_dir\n",
                        "cd \$INPUT_DIR\n",
                        "\n",
                        "\$BLAST_DIR/$type ",
                        "-query \$INPUT_DIR/split.$i/split.$i.fasta ",
                        "-db $database ",
                        "-outfmt $outputFormat ",
                        "-evalue $eval ",
                        "-num_threads 1 ",
                        "-max_target_seqs 1 ",
                        "\> $outputDir/split.$i.$type",
                        "\n",
                        "exit\n";
}
}
#print the task-array scripts hmmscan
if(! -e $hmmoutputDir){
        `mkdir $hmmoutputDir`;
}

for(my $i=1;$i<=$split_count;$i++){
     open my $HMM_SCRIPTS, ">", "$hmm-part-$i.sh" or die();
     print $HMM_SCRIPTS "#!/bin/bash\n",
                        "#PBS -N out.$hmm.part-$i\n",
                        "#PBS -l nodes=1:ppn=1,pmem=10G,pvmem=12g\n",
                        "#PBS -j oe\n",
                        "#PBS -l walltime=30:00:00\n",
                        "#PBS -M aadas\@uvm.edu\n",
                        "#PBS -m bea\n",
                        "\n",
                        "INPUT_DIR=$input_dir\n",
                        "HMM_DIR=$hmm_dir\n",
                        "cd \$INPUT_DIR\n",
                        "\n",
                        "\$HMM_DIR/$hmm ",
                        "--cpu 1 ",
                        "--domtblout $hmmoutputDir/split.$i.domtblout ",
                        "$hmmDB_dir/Pfam-A.hmm ",
                        "\$INPUT_DIR/split.$i/split.$i.fasta ",
                        "\n",
                        "exit\n";
}
#and submit it
```

Then (**When you run the 3rd script; generally it's submitting 300 jobs to VACC but it's not going to run; So, comment out blast-then it's only submitting hmm scan; Even if doesn't work then edit blast-part* and manually put part-1* then part-2*  **)

```
#!/bin/bash
cd `pwd`

perl batch_blastp_hmmscan.pl -q /users/a/a/aadas/nassellaBrachy_drought_freezing/blastP/npulBdis_Trinity211.fasta.transdecoder_dir/longest_orfs.pep -p blastp -d /users/a/a/aadas/blast_Brachyleytrum/database/uniprot_sprot.pep

chmod 700 *.sh
for i in blastp-part-2*; do
  qsub $i &
done

for i in hmmscan-part*; do
  qsub $i &
done
```

Run-

```
./run_blastp_hmmscan.sh
```

### Final step: predict the **likely coding regions**

```
#!/bin/bash
#PBS -N out.finalstep1
#PBS -l nodes=1:ppn=1,pmem=10G,pvmem=12g
#PBS -j oe
#PBS -l walltime=12:00:00
#PBS -M aadas@uvm.edu
#PBS -m bea

export PATH="/users/a/a/aadas/Bin/TransDecoder-3.0.1/transdecoder_plugins/cdhit:$PATH"
export PATH="/users/a/a/aadas/Bin/TransDecoder-3.0.1:$PATH"
export PATH="/users/a/a/aadas/Bin/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH"
export PATH="/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin:$PATH"

transDecoder_dir=/users/a/a/aadas/Bin/TransDecoder-3.0.1
INPUT_DIR=/users/a/a/aadas/blast_Brachyleytrum
cd $INPUT_DIR
######################################################################
###concatenate outputs for blastp and hmmscan searches
#####################################################################
cat blast_out/split.* > blastp.outfmt6
cat hmmscan_out/* > pfam.domtblout
####################################################################################################
####remove files that are no longer needed
###################################################################################################
rm -r split.* blastp-part-* hmmscan-part-*
#####################################################################################################
###submit final step of TransDecoder searching for potential coding regions of the transcripts
#####################################################################################################
$transDecoder_dir/TransDecoder.Predict -t $INPUT_DIR/Brachyleytrum_trinityv211.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
```

------

<div id='id-section9'/>

### Page 9: 2017-05-02. Principal Component Analysis (PCA)



------------------













​			
​		
​				
​		
​				
​		



