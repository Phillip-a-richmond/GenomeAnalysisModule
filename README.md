# Genome Analysis Module
> This course is intended for trainees in the Genetic Counselling MSc Program at the University of British Columbia. The computing resources are hosted at BC Children's Hospital Research Institute.  
> For more information contact Dr. Phillip Richmond (prichmond@bcchr.ca). If you are looking for content related to the manuscript, please see Release 1. 

---

## General Course Information & Prerequisites
### Course Format

> The course will be designed in a flexible format, adaptable for teaching online, or in person, or a combination of both. 

### Prerequisites  
Setting up the working environment (Video 1, 10 min). 
+ Content
  + Installing necessary software (see below) 
  + Connecting to the cluster 
  + Making a workspace 
+ Slides: [Link](https://docs.google.com/presentation/d/1tNzW21k7WnjfU-gDvL_gJB5R32GMlCrH0P1lHpH1jCY/edit#slide=id.g9d2975bcf0_0_5)

#### Software Installs  
+ For mac/linux users, only need native terminal which comes with the operating system.  You can access it in the Applications section for Mac users.  You also need to have X windows installed [Link](https://www.xquartz.org/)
+ For PC users, download and install [MobaXterm](http://mobaxterm.mobatek.net/) 
  + [Advanced MobaXterm usage](https://www.youtube.com/watch?v=Gkl8LD1rwlU)  
+ IGV installed on local machine 
  + [IGV install](https://www.broadinstitute.org/software/igv/log-in) 
+ Filezilla (or other file transfer client) installed on local machine
  + [FileZilla](https://filezilla-project.org/)
+ BCCHR VPN  


### Resources  
+ [Linux/Unix Cheatsheet](https://github.com/Phillip-a-richmond/ARC-Bioinformatics-Training/blob/master/Resources/UnixCheatSheet.pdf) 
+ [Linux/Unix online tutorial](http://www.ee.surrey.ac.uk/Teaching/Unix/) 
+ Editor cheat sheets 
  + [Emacs cheet sheet](http://www.rgrjr.com/emacs/emacs_cheat.html) 
  + [vi cheat sheet](http://www.lagmonster.org/docs/vi.html) 
  + [nano cheat sheet](http://www.codexpedia.com/text-editor/nano-text-editor-command-cheatsheet/) 
+ File transfer programs  
  + [FileZilla](https://filezilla-project.org/)
  + [WinSCP](https://winscp.net/eng/download.php) 
      
### To-Do before the first class  
Complete the Prerequisites (above). If you have any challenges, please contact Dr. Phillip Richmond (prichmond@bcchr.ca). 

---

## Course Syllabus  

### Prerequisite – to be done before Session 1 (15 minutes): 
Setting up the working environment 
- Installing secure shell and sshfs 
- Connecting to the cluster 
- Making a workspace
- [Link](https://docs.google.com/presentation/d/1tNzW21k7WnjfU-gDvL_gJB5R32GMlCrH0P1lHpH1jCY/edit#slide=id.g9d2975bcf0_0_5)
 
### Thursday June 3rd (10:00-12:00) - Intro to course, Intro to Genome Analysis for Rare Disease, Case Distribution 
- Soft-start (15 min: 10:00-10:15)
- Course Introduction (Lecture, 60 min: 10:15-11:15)
    - [Link](https://docs.google.com/presentation/d/1uzLQoHP3NODOWEs0BM2UZ8zG9lULiw25hcf07fLmSZs/edit#slide=id.p)
    - Goals of the course
    - Format of the course
    - Genome Analysis for Rare Disease Diagnosis
    - Limitations of whole exome sequencing, short reads, gene discovery, where the field is going
- Distribute cases & link to resources (Work-along, 30 min - 11:15-11:45)
    - Introduce concepts presented in the cases
    - OMIM, ClinVar, HPO
- Establish working environment (Work-along, 15 minutes - 11:45-12:00)
    - Install IGV
    - Login to cluster
    - Install SSHFS

### Friday June 4th (10:00-12:00) - Intro to linux command line & file handling
- Introduction to command line and HPC (Work-along, 45 min - 10:00-10:45)
    - [Link](https://docs.google.com/presentation/d/1TAIXbS3cJwEYpUXrhjig1e1JcVNYUvr8ZE0qWApuCA4/edit#slide=id.p)
    - Key commands for navigating
    - ls, mkdir, cd, pwd
    - Key commands for file handling (& exposure to genomics file formats)
    - more, less, grep, head, cut, tail, wc, chmod
    - Fasta, Fastq, VCF
- Introduction to file editing, bash scripts (Work-along, 30 min - 10:45-11:15)
    - Nano
    - Bash scripts
    - Variables
- *Problem Set 1 (45 min - 11:15-12:00) - basic linux navigation + file handling*


### Monday June 7th (10:00-12:00) - Short read mapping & visualization (Command line/bash)
- *Problem Set 1 Review (20 min: 10:00-10:20)*
- Short read mapping & Visualization (Work-along, 60 min - 10:20-11:20)
    - [Link](https://docs.google.com/presentation/d/1dDV7rpvDun7or578Pmm4fp9GRKRzMHKYpMs7AAlXZ0E/edit#slide=id.p)
    - Pipeline Overview
    - Exploring command line tools and how to use them
    - Map short reads to reference genome with BWA mem
    - Convert, sort, index, with Samtools
- *Problem Set 2 (40 min: 11:20-12:00) - short read mapping with bash scripts*


### Tuesday June 8th (10:00-12:00) - Process cases through pipeline using cluster
- *Problem Set 2 Review (30 min: 10:00-10:30)*
- Interacting with the cluster scheduler (Work-along, 45 min - 10:30-11:15)
    - Concept of job scheduler on cluster
    - Basic job scripts
    -  Mapping and visualization through cluster
- *Problem Set 3 (45 min - 11:15-12:00) -  Case processing fastq-->BAM + Visualize*

### Wednesday June 9th (10:00-12:00) - Variant calling and visualization
- *Problem Set 3 Review (30 min - 10:00-10:30)*
- Call variants with DeepTrio (Work-along, 45 min - 10:30-11:15) 
    - [Link](https://docs.google.com/presentation/d/1mLunTcJn1UIHSoHkBycS_6lWzvi_DGF9nNDSoo9tQ_g/edit#slide=id.p)
    - Call variants using pre-built DeepTrio scripts on pre-processed small BAM files
    - Visualize variants in IGV
    - Load GRCh38 genome
    - Load read files
    - Zoom to region with read coverage
    - Take snapshot of variant
- *Problem Set 4 (45 min - 11:15-12:00) - Call variants on samples*

### Thursday June 10th (10:00-12:00) - Variant annotation and interpretation
- *Extra time for Problem Set 4 (60 min - 10:00-11:00)*
- Variant interpretation for rare disease diagnosis (Lecture, 60 min - 11:00-12:00) 
    - Variant databases for annotation (ClinVar, gnomAD)
    - In silico metrics (CADD)
    - Phenotype:Genotype
    - Human Phenotype Ontology (HPO), 
    - Online Mendelian Inheritance in Man (OMIM)
    - Automation with Exomiser - Exomiser paper

### Friday June 11th (10:00-12:00) - Variant prioritization with Exomiser
- Running Exomiser on a VCF (Work-along, 45 min - 10:00-10:45)
    - [Link](https://docs.google.com/presentation/d/11LTBnsUh-gKxIMTfmSvuCX6R9YWsctNNMvLtUcu2Jhk/edit#slide=id.p)
    - Run exomiser using pre-built configs and scripts on Case 6
    - Interpreting Exomiser output in excel/html (Work-along, 45 min - 10:45-11:30)
    - Download exomiser output in TSV format, open and explore tabs in excel
    - Visualize output HTML file
- *Problem Set 5 (30 min - 11:30-12:00) - Run Exomiser on your case data*

### Monday June 14th (10:00-12:00) - Catch-up & open help session
- *Review Problem Set 5, address any issues, make sure everyone has their case exomiser output. (60 min - 10:00-11:00)*
- Open help session, working on Practical (60 min - 11:00-12:00)


### Tuesday June 15th (10:00-12:00) - Practical
- Practical (120 min - 10:00-12:00)
    - Each student presents their case to the group, 8 students ~15 minutes / presentation (10 min present, 5 min questions). 



