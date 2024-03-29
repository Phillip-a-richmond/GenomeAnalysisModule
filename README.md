# Genome Analysis Module
> This course is intended for trainees in the Genetic Counselling MSc Program at the University of British Columbia, and summer students within the Summer Student Research Program at BC Children's Hospital Research Institute. The computing resources are hosted at BC Children's Hospital Research Institute.  
> For more information contact Dr. Phillip Richmond (prichmond@bcchr.ca). If you are looking for content related to the manuscript, please see Release 1. 
> This version: November 2022


---

## General Course Information & Prerequisites
### Course Format

> The course will be designed in a flexible format, adaptable for teaching online, or in person, or a combination of both. 

### Prerequisites  
Setting up the working environment. 
+ Content
  + Installing necessary software (see below) 
  + Connecting to the cluster 
  + Making a workspace 
+ Slides: [Link](https://docs.google.com/presentation/d/10Nm2C2FyGBADgm2T_iKzF3sIB-wYWD1QMvgjLRCGXZQ/edit)

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
 
### Monday October 31st, 2022 (1:00-3:00) - Intro to course, Intro to Genome Analysis for Rare Disease, Case Distribution 
- Soft-start (15 min: 1:00-1:10)
- Course Introduction (Lecture, 75 min: 1:10-2:25)
    - [Link](https://docs.google.com/presentation/d/1JC4JiORk7zqOYeM-hK0vA-fmC2GdPNhY_UGbrOjMQFM/edit#slide=id.p)
    - Goals of the course
    - Format of the course
    - Genome Analysis for Rare Disease Diagnosis: [Link](https://docs.google.com/presentation/d/1nOmByd5QzS228ZbiuELQGQbiHnrAkSvkw4o9GNUEEfE/edit#slide=id.gdd01464b97_0_2279)
    - OMIM, ClinVar, HPO
    - Limitations of whole exome sequencing, short reads, gene discovery, where the field is going
- Distribute cases & link to resources (Work-along, 15 min - 2:25-2:40)
    - Introduce concepts presented in the cases
- Establish working environment (Work-along, 20 minutes - 2:40-3:00)
    - Install IGV
    - Login to cluster
    - Install SSHFS

### Wednesday November 2nd, 2022 (1:00-3:00) - Intro to linux command line & file handling
- Introduction to command line and HPC (Work-along, 45 min - 1:00-1:45)
    - [Link](https://docs.google.com/presentation/d/1zCOFq-lSbZu2oVnwTooATjWHhOfS3Vuw-aC1nou5ACQ/edit#slide=id.g17bc288be0f_0_411)
    - Key commands for navigating
    	- ls, mkdir, cd, pwd
    - Key commands for file handling (& exposure to genomics file formats)
    	- more, less, grep, head, cut, tail, wc, chmod
    - Fasta, Fastq, VCF
- Introduction to file editing, bash scripts (Work-along, 30 min - 1:45-2:15)
    - Nano
    - Bash scripts
    - Variables
- *Problem Set 1 (45 min - 2:15-3:00) - basic linux navigation + file handling*


### Friday November 4th, 2022 (8:30-10:30) - Short read mapping & visualization (Command line/bash)
- *Problem Set 1 Review (20 min: 8:30-8:50)*
- Short read mapping & Visualization (Work-along, 60 min - 8:50-9:50)
    - [Link](https://docs.google.com/presentation/d/1rdfC9FxQBya9EQpWPn1YhZ1XPvgNl-ek6Mrxauc-Kyc/edit)
    - Pipeline Overview
    - Exploring command line tools and how to use them
    - Map short reads to reference genome with BWA mem
    - Convert, sort, index, with Samtools
- *Problem Set 2 (40 min: 9:50-10:30) - short read mapping with bash scripts*


### Monday November 7th, 2022 (1:00-3:00) - Process cases through pipeline using cluster
- *Problem Set 2 Review (30 min: 1:00-1:30)*
- Interacting with the cluster scheduler (Work-along, 45 min - 1:30-2:15)
    - [Link](https://docs.google.com/presentation/d/1MQIV--4AwuI0Hl36HrzhR-buY_M8haX_zv5N78UJTNk/edit#slide=id.g165c9920962_0_132)
    - Concept of job scheduler on cluster
    - Basic job scripts
    -  Mapping and visualization through cluster
- *Problem Set 3 (45 min - 2:15-3:00) -  Case processing fastq-->BAM + Visualize*

### Wednesday November 9th, 2022 (1:00-3:00) - Variant calling and visualization
- *Problem Set 3 Review (30 min - 1:00-1:30)*
- Call variants with DeepTrio (Work-along, 45 min - 1:30-2:15) 
    - [Link](https://docs.google.com/presentation/d/11ou1MZ1WWrznPypNGyjB7SPKkw8mvNe9ny52IloWn7o/edit#slide=id.g10116e0bab2_0_35)
    - Call variants using pre-built DeepTrio scripts on pre-processed small BAM files
    - Visualize variants in IGV
    - Load GRCh38 genome
    - Load read files
    - Zoom to region with read coverage
    - Take snapshot of variant
- *Problem Set 4 (45 min - 2:15-3:00) - Call variants on samples*

### Thursday November 10th, 2022 (8:30-10:30) - Variant annotation and interpretation
- - *Extra time for Problem Set 4 (60 min - 8:30-9:30)*
- Open discussion/review of sessions 1-5
- Ensure that everyone has mapped, called variants, and visualized their case data


### Monday November 14th, 2022 (1:00-3:00) - Catch up and explain final activity
- [Link](https://docs.google.com/presentation/d/10KqxGE_l4lzCLmcHTNLpjEfUhPUXM5c3CtqUGPm5PgA/edit#slide=id.g165d2176f5a_0_132)
    - What is Exomiser?
    - Run exomiser using pre-built configs and scripts on Case 4
- Interpreting Exomiser output in excel/html (Work-along, 45 min - 1:45-2:30)
    - Download exomiser output in TSV format, open and explore tabs in excel
    - Visualize output HTML file
- *Problem Set 5 (30 min - 2:30-3:00) - Run Exomiser on case 3*


### Wednesday November 16th, 2022 (1:00-3:00) - Catch-up & open help session
- *Review Problem Set 5, address any issues, make sure everyone has their case exomiser output. (60 min - 1:00-2:00)*
- Open help session, working on Practical (60 min - 2:00-3:00)


### Friday November 18th, 2022 (8:30-10:30) - Practical
- Practical (120 min - 8:30-10:30)
    - Each student presents their case to the group, 8 students ~15 minutes / presentation (10 min present, 5 min questions). 



