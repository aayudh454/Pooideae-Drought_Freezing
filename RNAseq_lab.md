# RNA seq. Pipeline

### Author: Aayudh Das

Login info: **ssh aadas@bluemoon-user2.uvm.edu**

### Data structure (*Brachyleytrum aristosum* data)

**Treatment:** There are three replicates for "pre-treatment" (1-3x), three for 24 h cold shock (4-6y), and two for 6 wk cold (7-8z). Not the perfect data set but okay given the species.

## Table of contents    
* [Page 1: 2017-03-20](#id-section1). Moving files and trimming 

* [Page 2: 2017-03-23](#id-section2). Trimming Ba2-3x, Ba4-6y, Ba7-8z

* [Page 3 2017-03-24](#id-section3). Concatenation and assembly (using Trinity 2.4.0)

* [Page 4 2017-03-28](#id-section4). Assembly using Trinity 2.0.6 and 2.1.1 

* [Page 5 2017-04-11](#id-section5). Transcript quantification by RSEM

* [Page 6 2017-04-13](#id-section6). Build Transcript and Gene Expression Matrices

* [Page 7 2017-04-14](#id-section7). Differential Expression Analysis

* [Page 8 2017-04-25](#id-section8). Coding Region Identification in Trinity Assemblies

* [Page 9 2017-05-02](#id-section9). Principal Component Analysis (PCA)

  ​

------
<div id='id-section1'/>
###Page 1: 2017-03-20. Moving files, basics and trimming

#### 1. Where are our files?

After login you should **cd** space the location path. **ll** to see what are the files present

```
[aadas@bluemoon-user2 ~]$ cd /gpfs2/scratch/djshirle/MPS/170216_SNL128_0151_AHC72LBCXY/samples_out/
[aadas@bluemoon-user2 samples_out]$ ll
total 144
-rw-r--r-- 1 djshirle usr 2504 Feb 24 16:16 170216_demux_sheet-Preston.csv.quality_table_by_lane.txt
-rw-r--r-- 1 djshirle usr 2488 Feb 24 16:16 170216_demux_sheet-Preston.csv.quality_table_by_sample.txt
drwxr-xr-x 4 djshirle usr 8192 Feb 24 16:11 Ba1x
drwxr-xr-x 4 djshirle usr 8192 Feb 24 15:58 Ba2x
drwxr-xr-x 4 djshirle usr 8192 Feb 24 16:04 Ba3x
drwxr-xr-x 4 djshirle usr 8192 Feb 24 15:59 Ba4y
drwxr-xr-x 4 djshirle usr 8192 Feb 24 16:06 Ba5y
drwxr-xr-x 4 djshirle usr 8192 Feb 24 16:00 Ba6y
drwxr-xr-x 4 djshirle usr 8192 Feb 24 16:03 Ba7z
drwxr-xr-x 4 djshirle usr 8192 Feb 24 15:58 Ba8z
```

#### 2. How can you move all the files to your own directory from the server?

```
[aadas@bluemoon-user2 ~]$ cp -r /gpfs2/scratch/djshirle/MPS/170216_SNL128_0151_AHC72LBCXY/samples_out/* . &
```

cp=copy, *= means everything (for individual file just write the name of the file instead of *)

#### 3. You need to install trimmomatic in your PC

go to http://www.usadellab.org/cms/?page=trimmomatic

right click on binary and copy link address 

```
[aadas@bluemoon-user2 ~]$ wget -c http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
```

Now it's downloaded to your main directory as a zip file (Trimmomatic-0.36.zip)

you need to **unzip**

```
[aadas@bluemoon-user2 ~]$ unzip Trimmomatic-0.36.zip 
```

then you will see Trimmomatic-0.36 in your directory

#### 4. Create a new folder (command=mkdir) where you will execute all your analysis

```
[aadas@bluemoon-user2 ~]$ mkdir Ba
```

my new folder is Ba

#### 5. Let's see what files are present in Ba1x folder "pre-treatment" (1x)

```
[aadas@bluemoon-user2 ~]$ cd Ba1x/
[aadas@bluemoon-user2 Ba1x]$ ll
total 10661240
drwxr-xr-x 4 aadas usr        512 Mar 20 11:57 jcpresto_BrArisVR_20170217_Ba1x_R1_fastqc
-rw-r--r-- 1 aadas usr     197928 Mar 20 11:57 jcpresto_BrArisVR_20170217_Ba1x_R1_fastqc.zip
-rw-r--r-- 1 aadas usr 5399685374 Mar 20 11:56 jcpresto_BrArisVR_20170217_Ba1x_R1.fastq.gz
-rw-r--r-- 1 aadas usr        157 Mar 20 11:56 jcpresto_BrArisVR_20170217_Ba1x_R1.fastq.gz.md5sum
drwxr-xr-x 4 aadas usr        512 Mar 20 11:57 jcpresto_BrArisVR_20170217_Ba1x_R2_fastqc
-rw-r--r-- 1 aadas usr     201979 Mar 20 11:57 jcpresto_BrArisVR_20170217_Ba1x_R2_fastqc.zip
-rw-r--r-- 1 aadas usr 5517001528 Mar 20 11:57 jcpresto_BrArisVR_20170217_Ba1x_R2.fastq.gz
-rw-r--r-- 1 aadas usr        157 Mar 20 11:57 jcpresto_BrArisVR_20170217_Ba1x_R2.fastq.gz.md5sum
```

Now you need **R1.fastq.gz** and  **R2.fastq.gz** files, so copy them to your Ba folder. **../Ba/** = this means your are moving to Ba folder. **Ba1x_precold.R1.fastq.gz** is the new name.

```
[aadas@bluemoon-user2 Ba1x]$ cp jcpresto_BrArisVR_20170217_Ba1x_R1.fastq.gz ../Ba/Ba1x_precold.R1.fastq.gz
[aadas@bluemoon-user2 Ba1x]$ cp jcpresto_BrArisVR_20170217_Ba1x_R2.fastq.gz ../Ba/Ba1x_precold.R2.fastq.gz
```

Now in your Ba folder you have both the R1 and R2 files-

```
[aadas@bluemoon-user2 Ba]$ ll
total 10660840
-rw-r--r-- 1 aadas usr 5399685374 Mar 20 12:21 Ba1x_precold.R1.fastq.gz
-rw-r--r-- 1 aadas usr 5517001528 Mar 20 12:25 Ba1x_precold.R2.fastq.gz
```

#### 6. Now create a script for running the trimming program

you should be in your **Ba** folder 

```
[aadas@bluemoon-user2 ~]$ cd Ba/
[aadas@bluemoon-user2 Ba]$
```

 Now to create a blank script

```
[aadas@bluemoon-user2 ~]$ vi trimmomatic.sh
```

Now press **i** to insert 

Copy and paste everything present in the new script from this sample one

```
#!/bin/bash

######## This job needs 1 nodes, 2 processors total
#PBS -l nodes=1:ppn=2
# it needs to run for 6 hours
#PBS -l walltviime=06:00:00
#PBS -N renamer
#PBS -j oe
#PBS -M YOUR_ACCOUNT@uvm.edu
#PBS -m bea
###LOAD JAVA MODULE AVAILABLE FROM THE CLUSTER, YOU MAY WANT TO CHECK FIRST
module load java-sdk/sun-jdk-1.6.0.12
ulimit -s unlimited
###CHANGE THE DIRECTORY ACCORDINGLY, THE FOLLOWING SETTINGS ARE FOR MY ACCOUNT
SOFTWARE=/users/j/z/jzhong2/bin/Trimmomatic-0.33/
workDIR=/users/j/z/jzhong2/scratch/wheat_polyploidization
cd $workDIR
#####TRIMMING COMMANDS AND PARAMETERS
java -jar $SOFTWARE/trimmomatic-0.33.jar PE -phred33 $workDIR/Melica6weekcold01.R1.fq $workDIR/Melica6weekcold01.R2.fq $workDIR/Melica6weekcold01_R1.trimmo.fq $workDIR/Melica6weekcold01_R1.unpaired.fq $workDIR/Melica6weekcold01_R2.trimmo.fq $workDIR/Melica6weekcold01_R2.unpaired.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:40
```

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
#PBS -N renamer
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea
###LOAD JAVA MODULE AVAILABLE FROM THE CLUSTER, YOU MAY WANT TO CHECK FIRST
ulimit -s unlimited
###CHANGE THE DIRECTORY ACCORDINGLY, THE FOLLOWING SETTINGS ARE FOR MY ACCOUNT
SOFTWARE=/users/a/a/aadas/Trimmomatic-0.36
workDIR=/users/a/a/aadas/Ba
cd $workDIR
#####TRIMMING COMMANDS AND PARAMETERS
java -jar $SOFTWARE/trimmomatic-0.36.jar PE -phred33 $workDIR/Ba1x_precold.R1.fastq.gz $workDIR/Ba1x_precold.R2.fastq.gz $workDIR/Ba1x_precold.R1.trimmo.fq.gz $workDIR/Ba1x_precold.R1.unpaired.fq.gz $workDIR/Ba1x_precold.R2.trimmo.fq.gz $workDIR/Ba1x_precold.R2.unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:40
```

i. You need to check **Submitting Jobs to the Cluster**-

https://www.uvm.edu/~vacc/?Page=userguide.php

first part of the script explains that details

ii. Now you need to specify where your software is present i.e. the Trimmomatic-0.36 which is in your main directory **/users/a/a/aadas/Trimmomatic-0.36**

iii. Specify your working directory-**/users/a/a/aadas/Ba** (because your R1 and R2 files are in Ba)

iv. TRIMMING COMMANDS AND PARAMETERS

a. change software version from as **trimmomatic-0.36**

b. Now the **first 2 files are your input file**, so after $workDIR/"name of the file" space [here R1 and R2 is the main change]

c. Last 4 files are **output files**.  

$workDIR/Ba1x_precold_R1.trimmo.fq.gz $workDIR/Ba1x_precold.R1.unpaired.fq.gz $workDIR/Ba1x_precold.R2.trimmo.fq.gz $workDIR/Ba1x_precold.R2.unpaired.fq.gz

for 1st pair one is paired named as R1.**trimmo** and R1.**unpaired**. same for the second pair.

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
### Page 2: 2017-03-23. Trimming Ba2-3x, Ba4-6y, Ba7-8z 

###### Some moving tips 

#### Rename a file

```
mv oldname newname
mv *.trimmo.fq.gz ../Brachyleytrum_aristosum
```

#### Moving Ba2x,Ba3x file to Ba folder

```
[aadas@bluemoon-user2 ~]$ cd Ba2x/
[aadas@bluemoon-user2 Ba2x]$ cp jcpresto_BrArisVR_20170217_Ba2x_R1.fastq.gz ../Ba/Ba2x_precold.R1.fastq.gz
[aadas@bluemoon-user2 Ba2x]$ cp jcpresto_BrArisVR_20170217_Ba2x_R2.fastq.gz ../Ba/Ba2x_precold.R2.fastq.gz

```

#### Now edit the script -

```
[aadas@bluemoon-user2 ~]$ cd Ba/
[aadas@bluemoon-user2 Ba]$ ll
[aadas@bluemoon-user2 Ba]$ vi trimmomatic.sh 
```

replace the Ba1x with Ba2x. Then :wq to save

#### Now make the script execute the task 

```
[aadas@bluemoon-user2 Ba]$ chmod 700 trimmomatic.sh 
```

#### Submit the job and view

```
[aadas@bluemoon-user2 Ba]$ qsub trimmomatic.sh
[aadas@bluemoon-user2 Ba]$ showq -u aadas
```

#### Copy only the trimmo.fq.gz files to a new folder for assembly

```
[aadas@bluemoon-user2 Ba]$ cp *.trimmo.fq.gz ../Brachyleytrum_aristosum
```



------

<div id='id-section3'/>

### Page 3: 2017-03-24. Concatenation and assembly

###### Some tips

to delete a line: **esc** and **d+d**

https://www.cs.colostate.edu/helpdocs/vi.html

Save all the commands as a history.txt

```
[aadas@bluemoon-user2 ~]$ history > history.txt
```

#### 1. Install trinity (using Trinity 2.4.0) 

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

#### 2. Execute concatenation for both R1 and R2

```
[aadas@bluemoon-user2 Brachyleytrum_aristosum]$ zcat *R1.trimmo.fq.gz > BrachyletrumARI.R1.trimmo.fq &
[aadas@bluemoon-user2 Brachyleytrum_aristosum]$ zcat *R2.trimmo.fq.gz > BrachyletrumARI.R2.trimmo.fq &
```

#### After finishing the concatenation check the "sequence header for both R1 and R2"-it should be same

```
[aadas@bluemoon-user2 Brachyleytrum_aristosum]$ grep -c "@" BrachyletrumARI.R1.trimmo.fq 
[aadas@bluemoon-user2 Brachyleytrum_aristosum]$ grep -c "@" BrachyletrumARI.R2.trimmo.fq 
```

both shows 168097158; That means R1 and R2 has same reads.

#### 3. Make sure you have a script to do the assembly 

```
#!/bin/bash

#PBS -l nodes=1:ppn=8,mem=96G,vmem=100G
#PBS -q poolmemq
# it needs to run for 6 hours
#PBS -l walltime=30:00:00
#PBS -N trinity
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea
module load samtools-1.3.1-gcc-6.3.0-e5jw5u4
module load bowtie2-2.2.5-gcc-6.3.0-daskah5
ulimit -s unlimited

SOFTWAREDIR=/users/a/a/aadas/Bin/trinityrnaseq-Trinity-v2.4.0
WORKINGDIR=/users/a/a/aadas/Brachyleytrum_aristosum
cd $WORKINGDIR

$SOFTWAREDIR/Trinity --seqType fq --max_memory 96G --left $WORKINGDIR/AegilopsLON_reads1.fastq --right $WORKINGDIR/AegilopsLON_reads2.fastq --CPU 8
```

edit the script (change WORKINGDIR file names)

```
#!/bin/bash

#PBS -l nodes=1:ppn=16,mem=96G,vmem=100G
#PBS -q poolmemq
# it needs to run for 6 hours
#PBS -l walltime=30:00:00
#PBS -N trinity
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea
module load samtools-1.3.1-gcc-6.3.0-e5jw5u4
module load bowtie2-2.2.5-gcc-6.3.0-daskah5
ulimit -s unlimited

SOFTWAREDIR=/users/a/a/aadas/Bin/trinityrnaseq-Trinity-v2.4.0
WORKINGDIR=/users/a/a/aadas/Brachyleytrum_aristosum
cd $WORKINGDIR

$SOFTWAREDIR/Trinity --seqType fq --max_memory 96G --left $WORKINGDIR/BrachyletrumARI.R1.trimmo.fq --right $WORKINGDIR/BrachyletrumARI.R2.trimmo.fq --CPU 16
```

Here you need to include the module load but first you need to check what modules are availble in the cluster.

```
[aadas@bluemoon-user2 Brachyleytrum_aristosum]$ module avail
```

Then include the samtools and bowtie version, so that your script can use that program.

Increase the cpu load to 16 in **ppn=16** and also in the end **CPU 16** 

Submit the job and view

```
[aadas@bluemoon-user2 Brachyleytrum_aristosum]$ qsub vacctrinity.sh 
[aadas@bluemoon-user2 Brachyleytrum_aristosum]$ showq -u aadas
```

If you see that 30% completed after 30hrs then qsub again and it will catch up from where it was left before.

------

<div id='id-section4'/>

### Page 4: 2017-03-28. Assembly using Trinity 2.0.6

Check the assembly file

```
# number of sequence present
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta | less
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | less
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | sort -u | less
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | sort -u | wc -l
```

94102

check no. of seq

```
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta -c
```

567758

#### New script for assembly (Trinity 2.0.6) 

```
#!/bin/bash

#PBS -l nodes=1:ppn=16,mem=96G,vmem=128G
#PBS -q poolmemq
# it needs to run for 6 hours
#PBS -l walltime=30:00:00
#PBS -N trinity
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea
module load samtools-1.3.1-gcc-6.3.0-e5jw5u4

export PATH="/users/a/a/aadas/Bin/bowtie-1.2:$PATH"
export PATH="/users/a/a/aadas/Bin/jre1.7.0_80/bin:$PATH"
export PATH="/users/a/a/aadas/Bin/jre1.7.0_80/bin/java:$PATH"

ulimit -s unlimited

SOFTWAREDIR=/users/a/a/aadas/Bin_Trinity206/trinityrnaseq-2.0.6
WORKINGDIR=/users/a/a/aadas/Brachyleytrum_aristosum
cd $WORKINGDIR

$SOFTWAREDIR/Trinity --seqType fq --max_memory 96G --left $WORKINGDIR/BrachyletrumARI.R1.trimmo.fq --right $WORKINGDIR/BrachyletrumARI.R2.trimmo.fq --CPU 16
```

##### When you download something in mac and upload that to server

```
scp bowtie-1.2-linux-x86_64.zip aadas@bluemoon-user2.uvm.edu:~/
```

##### Delete an entire directory

```
rm -rf filename
```

### Running assembly with Trinity 2.1.1

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
WORKINGDIR=/users/a/a/aadas/Brachyleytrum_aristosum/Trinity211assembly
cd $WORKINGDIR

/users/a/a/aadas/Bin/trinityrnaseq-2.1.1/Trinity --seqType fq --normalize_reads --max_memory 256G --left /users/a/a/aadas/Brachyleytrum_aristosum/Trinity211assembly/BrachyletrumARI.R1.trimmo.fq --right /users/a/a/aadas/Brachyleytrum_aristosum/Trinity211assembly/BrachyletrumARI.R2.trimmo.fq --CPU 24
```

 

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

Check back after your job is done

```
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta -c
303494
[aadas@bluemoon-user2 trinity_out_dir]$ grep ">" Trinity.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | sort -u | wc -l
173604
```

l

------

<div id='id-section5'/>

### Page 5: 2017-04-11. Transcript quantification by RSEM

RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome

Script for transcript quantification.

* Make sure your Trinity.fasta file in the directory
* Make sure that R1.trimmo.fq.gz and R2.trimmo.fq.gz are also in the FASTQ_DIR

### For pre-cold

```
#!/bin/bash

######## This job needs 1 nodes, 4 processors total
#PBS -l nodes=1:ppn=4,pmem=8gb,pvmem=9gb
#PBS -l walltime=30:00:00
#PBS -N outprecold
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea
module load samtools-1.3.1-gcc-6.3.0-e5jw5u4

TRINITY_HOME=/users/a/a/aadas/Bin/trinityrnaseq-2.1.1

export PATH=/users/a/a/aadas/Bin/bowtie-1.1.1:$PATH
export PATH=/users/a/a/aadas/Bin/RSEM-1.2.19:$PATH

FASTQ_DIR=/users/a/a/aadas/Brachyleytrum_aristosum/assemblyTrinity2.1.1
cd $FASTQ_DIR

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts Brachyleytrum_trinityv211.fasta --seqType fq --left $FASTQ_DIR/Ba1x_precold_R1.trimmo.fq.gz --right $FASTQ_DIR/Ba1x_precold_R2.trimmo.fq.gz --est_method RSEM --aln_method bowtie --thread_count 4 --SS_lib_type RF --trinity_mode --prep_reference --output_dir $FASTQ_DIR --output_prefix Brachyleytrum_precold01


$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts Brachyleytrum_trinityv211.fasta --seqType fq --left $FASTQ_DIR/Ba2x_precold_R1.trimmo.fq.gz --right $FASTQ_DIR/Ba2x_precold_R2.trimmo.fq.gz --est_method RSEM --aln_method bowtie --thread_count 4 --SS_lib_type RF --trinity_mode --output_dir $FASTQ_DIR --output_prefix Brachyleytrum_precold02


$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts Brachyleytrum_trinityv211.fasta --seqType fq --left $FASTQ_DIR/Ba3x_precold_R1.trimmo.fq.gz --right $FASTQ_DIR/Ba3x_precold_R2.trimmo.fq.gz --est_method RSEM --aln_method bowtie --thread_count 4 --SS_lib_type RF --trinity_mode --output_dir $FASTQ_DIR --output_prefix Brachyleytrum_precold03
```

The above `R` commands list the top **10 highest expressed genes** (shown below).

```
data = read.table("Brachyleytrum_precold01.genes.results", header=T, stringsAsFactors=F)
idx = order(data[,"TPM"], decreasing=T)
data[idx[1:10], c("gene_id", "expected_count", "TPM")]
```

​                 gene_id expected_count      TPM

83500   TRINITY_DN32640_c2_g2      195870.99 88722.96
135779   TRINITY_DN5885_c0_g1       18720.00 21437.99
88514   TRINITY_DN33238_c0_g3        2823.75 15100.72
48900   TRINITY_DN25961_c0_g1       24537.02 14854.18
85085   TRINITY_DN32831_c0_g2       36801.00 14187.48
135853  TRINITY_DN58913_c0_g1       12804.00 13984.62
83499   TRINITY_DN32640_c2_g1       28301.73 12804.37
89635  TRINITY_DN33369_c0_g15       24298.99 11610.17
88512   TRINITY_DN33238_c0_g1        2153.17 11514.64
14964   TRINITY_DN16119_c1_g1       21600.00 10408.22

###  

------

<div id='id-section6'/>

### Page 6: 2017-04-13. Build Transcript and Gene Expression Matrices

Terms:

**<u>FPKM</u>: fragments per kilobase transcript length per million fragments mapped**

**<u>TPM</u>: transcripts per million transcripts**

Using the transcript and gene-level abundance estimates for each of your samples, construct a matrix of counts and a matrix of normalized expression values using the following script:

#### For genes

```
[aadas@bluemoon-user2 assemblyTrinity2.1.1]$ ~/Bin/trinityrnaseq-2.1.1/util/abundance_estimates_to_matrix.pl --est_method RSEM Brachyleytrum_precold01.genes.results Brachyleytrum_precold02.genes.results Brachyleytrum_precold03.genes.results Brachyleytrum_coldshock01.genes.results Brachyleytrum_coldshock02.genes.results Brachyleytrum_coldshock03.genes.results Brachyleytrum_sixweekcold01.genes.results Brachyleytrum_sixweekcold02.genes.results --out_prefix Brachyleytrum.genes 
```

#### For isoforms

```
~/Bin/trinityrnaseq-2.1.1/util/abundance_estimates_to_matrix.pl --est_method RSEM Brachyleytrum_precold01.isoforms.results Brachyleytrum_precold02.isoforms.results Brachyleytrum_precold03.isoforms.results Brachyleytrum_coldshock01.isoforms.results Brachyleytrum_coldshock02.isoforms.results Brachyleytrum_coldshock03.isoforms.results Brachyleytrum_sixweekcold01.isoforms.results Brachyleytrum_sixweekcold02.isoforms.results --out_prefix Brachyleytrum.isoforms
```

#### Counting Numbers of Expressed Transcripts or Genes

```
~/Bin/trinityrnaseq-2.1.1/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \ Brachyleytrum.genes.TPM.not_cross_norm | tee Brachyleytrum.genes.TPM.not_cross_norm.counts_by_min_TPM

~/Bin/trinityrnaseq-2.1.1/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
          trans_matrix.TPM.not_cross_norm | tee trans_matrix.TPM.not_cross_norm.counts_by_min_TPM
```

The above table indicates that we have 52,624 'genes' that are expressed by at least 1 TPM in any one of the many samples in this expression matrix.

Plotting the number of 'genes' (or 'transcripts') as a function of minimum TPM threshold, we can see that the vast majority of all expressed features have very little expression support. Using R (or your own favorite data analysis package), we might extrapolate the number of expressed 'genes' based on the trend prior to the massive influx of lowly expressed transcripts:

```
setwd("~/Desktop")
list.files()
data = read.table("Brachyleytrum.genes.TPM.not_cross_norm.counts_by_min_TPM", header=T)
plot(data, xlim=c(-100,0), ylim=c(0,100000), t='b')
```

![alt tag](https://github.com/aayudh454/Lab-RNA-seq-pipeline-/blob/master/Rplot.tiff)

![Rplot](/Users/aayudhdas/Desktop/Lab-RNA-seq-pipeline-/Rplot.tiff)

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

![Rplot01](/Users/aayudhdas/Desktop/Lab-RNA-seq-pipeline-/Rplot01.tiff)

The linear regression allows us to extrapolate (based on the Y-intercept) that we have 19361 'genes', which is a far better guess than our count of 52,624 'genes' having at least 1 TPM in any sample, and certainly better than the 1.4 million 'genes' that were assembled. 

------

<div id='id-section7'/>

### Page 7: 2017-04-14. Differential Expression Analysis

First lets get 'R' working

```
q[aadas@bluemoon-user2 ~]$ module load r-3.3.2-gcc-6.3.0-bmdvb4s
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

### Identifying DE Features: No Biological Replicates 	

### First lets get 'R' working		

```
module load r-3.3.2-gcc-6.3.0-bmdvb4s
```

It's very important to have biological replicates to power DE detection and reduce false positive predictions. If you do not have biological replicates, edgeR will allow you to perform DE analysis if you manually set the --dispersion parameter. Values for the **dispersion parameter** must be chosen carefully, and you might begin by exploring values between **0.1 and 0.4**.

now, run edgeR via the helper script provided in the Trinity distribution:

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Brachyleytrum.genes.counts.matrix --method edgeR --dispersion 0.1 --output Brachyleytrum_edgeR
```

**Brachyleytrum_coldshock01_vs_Brachyleytrum_coldshock02**

![alt tag] (https://github.com/aayudh454/Lab-RNA-seq-pipeline-/blob/master/Brachyleytrum.genes.counts.matrix.Brachyleytrum_coldshock01_vs_Brachyleytrum_coldshock02.edgeR.DE_results.MA_n_Volcano.pdf)

### Identifying DE features: With biological replicates (PREFERRED)

Be sure to have a tab-delimited 'samples_described.txt' file that describes the relationship between samples and replicates. For example:

just create a text file and copy paste this

```
conditionA   Brachyleytrum_precold01
conditionA   Brachyleytrum_precold02
conditionA   Brachyleytrum_precold03

conditionB   Brachyleytrum_coldshock01
conditionB   Brachyleytrum_coldshock02
conditionB   Brachyleytrum_coldshock03

conditionC   Brachyleytrum_sixweekcold01
conditionC   Brachyleytrum_sixweekcold02
```

#### voom method

Any of the available methods support analyses containing biological replicates. Here, for example, we again choose voom within the limma package.

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Brachyleytrum.genes.counts.matrix --method voom --samples_file samples_described.txt  
```

![alt tag](https://github.com/aayudh454/Lab-RNA-seq-pipeline-/blob/master/Brachyleytrum.genes.counts.matrix.conditionA_vs_conditionB.voom.DE_results.MA_n_Volcano.pdf)

### Interactive Volcano and MA Plots using Glimma

The [Glimma](https://bioconductor.org/packages/release/bioc/html/Glimma.html) software provides interactive plots. Generate volcano and MA-plots for any of your pairwise DE analysis results like so:

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/glimma1.R --samples_file samples_described.txt --DE_results Brachyleytrum.genes.counts.matrix.conditionA_vs_conditionB.voom.DE_results --counts_matrix Brachyleytrum.genes.counts.matrix
```

### didn't work!

### Extracting and clustering differentially expressed transcripts

An initial step in analyzing differential expression is to extract those transcripts that are most differentially expressed (most significant FDR and fold-changes) and to cluster the transcripts according to their patterns of differential expression across the samples. 

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix Brachyleytrum.genes.counts.matrix -P 1e-3 -C 2 --samples samples_described.txt
```

you might need to update the q value in library.

which will extract all genes that have P-values at most 1e-3 and are at least 2^2 fold differentially expressed. For each of the earlier pairwise DE comparisons, this step will generate the following files:

## Gene Ontology (GO) Enrichment Analysis on Differentially Expressed Genes

There are three different methods for partitioning genes into clusters:

- use K-means clustering to define K gene sets. (use the -K parameter). This does not leverage the already hierarchically clustered genes as shown in the heatmap, and instead uses a least-sum-of-squares method to define exactly k gene clusters.
- cut the hierarchically clustered genes (as shown in the heatmap) into exactly K clusters.
- (Recommended) cut the hierarchically clustered gene tree at --Ptree percent height of the tree.

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -R diffExpr.P1e-3_C2.matrix.RData --Ptree 60
```

A directory will be created called: 'diffExpr.P0.001_C2.matrix.RData.clusters_fixed_P_60' that contains the expression matrix for each of the clusters (log2-transformed, median centered). A summary pdf file is provided as 'my_cluster_plots.pdf' that shows the expression patterns for the genes in each cluster. Each gene is plotted (gray) in addition to the mean expression profile for that cluster (blue), as shown below:



------

<div id='id-section8'/>

### Page 8: 2017-04-25. Coding Region Identification in Trinity Assemblies

Follow this link-http://transdecoder.github.io/

The *TransDecoder* utility is run on a fasta file containing the target transcript sequences. The simplest usage is as follows:

Step 1: extract the **long open reading frames**

```
[aadas@bluemoon-user2 annotation]$ /users/a/a/aadas/Bin/TransDecoder-3.0.1/TransDecoder.LongOrfs -t Brachyleytrum_trinityv211.fasta
```

Step 2 predict the **likely coding regions**

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

### Step 1

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
INPUT_DIR=/users/a/a/aadas/blast_Brachyleytrum
cd $INPUT_DIR

$transDecoder_dir/TransDecoder.LongOrfs -t $INPUT_DIR/Brachyleytrum_trinityv211.fasta
```

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
my $hmmDB_dir="/users/a/a/aadas/blast_Brachyleytrum/database";
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

perl batch_blastp_hmmscan.pl -q Brachyleytrum_trinityv211.fasta.transdecoder_dir/longest_orfs.pep -p blastp -d /users/a/a/aadas/blast_Brachyleytrum/database/uniprot_sprot.pep

chmod 700 *.sh
for i in blastp-part*; do
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

### Final step

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

Another important analysis method to explore relationships among the sample replicates is Principal Component Analysis (PCA). You can generate a PCA plot like so:

```
[aadas@bluemoon-user2 pcaplots]$ module load r-3.3.2-gcc-6.3.0-bmdvb4s
```

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/PtR --matrix Trinity_trans.counts.matrix \
    -s samples.txt --log2 --prin_comp 3
```

```
 %  ~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/PtR --matrix Brachyleytrum.genes.counts.matrix \
    -s samples_described.txt --log2 --prin_comp 3
```

The --prin_comp 3 indicates that the first three principal components will be plotted, as shown above, with PC1 vs. PC2 and PC2 vs. PC3. In this example, the replicates cluster tightly according to sample type, which is very reassuring.

If you have replicates that are clear outliers, you might consider removing them from your study as potential confounders. If it's clear that you have a [batch effect](http://www.nature.com/nrg/journal/v11/n10/full/nrg2825.html), you'll want to eliminate the batch effect during your downstream analysis of differential expression.

------------------

<div id='id-section10'/>

### Page 10: 2017-10-02. Orthofinder setup

https://github.com/davidemms/OrthoFinder

#### some tips (If you want to upload a fle from desktop)

first cd to Desktop

```
AAYUDHs-MacBook-Air:Desktop aayudhdas$ scp fastme-2.1.5.tar.gz aadas@bluemoon-user2.uvm.edu:~/Bin
```

#### 1. Get **MCL-edge software **

https://micans.org/mcl/

After you tar -zxvf

```
cd mcl-12-068
./configure --prefix=$HOME/local
make install
```

#### 2. FastME 2.0

```
./configure
make
```



#### 3. DLCpar

https://www.cs.hmc.edu/~yjw/software/dlcpar/



-------

### Extra

Follow this link-http://trinotate.github.io/

Trinotate **relies heavily on SwissProt and Pfam**, and custom protein files are generated as described below to be specifically used with Trinotate. You can obtain the protein database files by running this Trinotate build process. This step will download several data resources including the latest version of swissprot, pfam, and other companion resources, create and populate a Trinotate boilerplate sqlite database (Trinotate.sqlite), and yield *uniprot_sprot.pep* file to be used with BLAST, and the *Pfam-A.hmm.gz* file to be used for Pfam searches. Run the build process like so:

```
[aadas@bluemoon-user2 Trinotate-3.0.2]$ ~/Bin/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
```

and once it completes, it will provide to you:

```
Trinotate.sqlite
uniprot_sprot.pep
Pfam-A.hmm.gz
```

Prepare the protein database for blast searches by:	

```
[aadas@bluemoon-user2 Trinotate-3.0.2]$ ~/Bin/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl makeblastdb -in uniprot_sprot.pep -dbtype prot
```

script to do the analysis

```
#!/bin/bash

######## This job needs 1 nodes, 4 processors total
#PBS -l nodes=1:ppn=4,pmem=8gb,pvmem=9gb
#PBS -l walltime=30:00:00
#PBS -N blastp
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea

# This single line using the blastp command below will compare your transcript fasta file
# (-query) to the already formatted uniref90 database (-db).
# You can enter 'blastp --help' for a list of the parameters.
# We choose the tab-delimited output format (6) and to only help the top hit (-max_target_seqs)
# and only if it has a minimum evalue of 0.001.

export PATH="/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin:$PATH"

blastp -query /users/a/a/aadas/annotation/blastp/Brachyleytrum_orfs.pep \
       -db /users/a/a/aadas/annotation/blastp/uniprot_sprot.pep \
       -out /users/a/a/aadas/annotation/blastp/blastp_vs_uniprot.outfmt6 \
       -outfmt 6 \
       -evalue 1e-3 \
       -max_target_seqs 1
```















​			
​		
​				
​		
​				
​		



