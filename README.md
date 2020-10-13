# Genome Analysis Module
> This course is intended for trainees in the Genetic Counselling MSc Program at the University of British Columbia. The computing resources are hosted at BC Children's Hospital Research Institute.  
> For more information contact Dr. Phillip Richmond (prichmond@bcchr.ca). If you are looking for content related to the manuscript, please see Release 1. 

---
## Manuscript
> A manuscript detailing the prior work can be found here:

[Introduction to Genomic Analysis Workshop: A catalyst for engaging life-science researchers in high throughput analysis](https://f1000research.com/articles/8-1221)


## General Course Information & Prerequisites
### Course Format

> The course will be designed in a flexible format, adaptable for teaching online, or in person, or a combination of both. 

### Prerequisites  
Setting up the working environment (Video 1, 10 min). 
+ Content
  + Installing necessary software (see below) 
  + Connecting to the cluster 
  + Making a workspace 
+ Slides: 
+ Video: 
+ Worksheet: 

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

### Session 1: Monday, November 2nd, 1:00PM-5:00PM (PST)
+ Course Introduction (Lecture: 1:00-1:45)
  + Goals of the course
  + Format of the course
  + Genome Analysis for Rare Disease Diagnosis 
  + Slides[]()
+ Introduction to command line and HPC (Video 2, 15-20 min: 1:45-2:45)
  + Key commands for navigating: cp, ls, mv, cut, mkdir, rm, scp, wget
  + Key commands for file handling (& exposure to genomics file formats): more, less, grep, head, cut, tail, wc
  + Key file formats: Fasta, Fastq, VCF
  + Slides:[]()
  + Video:[]()
+ Coffee Break and catch up (2:45-3:00)
+ Introduction to file editing, bash scripts (Video 3, 15-20 min: 3:00-4:00)
  + Emacs
  + Bash scripts
  + Slides:[]()
  + Video:[]()
+ *Problem Set 1* (4:00-5:00)
  + Link: []()

### Session 2: Wednesday, November 4th, 9:00AM-12:00PM (PST)
+ *Problem Set 1 Review* (9:00-10:00)
+ Interacting with the cluster scheduler (Lecture: 10:00-10:30)
  + Learning the cluster
  + SLURM job submitting
  + Module system for software
  + Slides: []()
  + Video: []()
+ Interacting with the cluster scheduler (Video 4, 15 min: 10:30-11:00)
  + Edit job scripts, submit to queue, collect output
+ Coffee Break and catch up (11:00-11:15)
+ *Problem Set 2* (11:15-12:00)
  + Link: []()

### Session 3: Wednesday, November 4th, 2:00PM-3:45PM (PST)
+ *Problem Set 2 Review* (2:00-2:30)
+ Short read mapping and visualization (Lecture: 2:30-2:45)
  + Pipeline Overview
  + Exploring command line tools and how to use them
+ Mapping Reads to the genome (Video 5, 15 min: 2:45-3:15) 
  + Map short reads with BWA mem
  + Convert SAM to BAM using Samtools
  + Use Samtools flagstat to get read mapping statistics
+ *Problem Set 3* (3:15-3:45)
  + Link: []() 

### Session 4: Monday, November 9th, 1:00PM-5:00PM (PST)
+ *Problem Set 3 Review* (1:00-2:00)
+ Variant calling and visualization (Lecture: 2:00-2:30)
  + Many tools for variant calling, but require same inputs
  + Different classes of variants which can be called 
  + Call variants with DeepVariant (Video 6, 15 min: 2:30-3:00) 
  + Call variants using pre-built DeepVariant scripts on pre-processed BAM files
  + Slides: []()
+ Coffee Break and catch up (3:00-3:15)
+ Visualization with IGV (Video 7, 10 min: 3:15-4:00)
  + Download and install IGV
  + Load GRCh38 genome
  + Load read files
  + Zoom to region with read coverage
  + Take snapshot of variant
  + Slides: []()
  + Video: []()
+ *Problem Set 4* (4:00-5:00)

_*Any assistance needed with Problem Sets 1-4 covered in the week between sessions. At this point, we need to ensure that all students have processed their samples to the VCF stage.*_

### Session 5: Monday, November 16th, 1:00PM-5:00PM (PST), 
+ *Problem Set 4 Review*
+ Variant interpretation for rare disease diagnosis (Lecture: 2:00-2:45) 
  + Variant databases for annotation (ClinVar, gnomAD)
  + In silico metrics (CADD)
  + Phenotype:Genotype
    + Human Phenotype Ontology (HPO), 
    + Online Mendelian Inheritance in Man (OMIM)
  + Automation with Exomiser
  + Slides: []()
+ Running Exomiser on a VCF (Video 9, 15 min: 2:45-3:15)
  + Run exomiser using pre-built configs and scripts on pre-organized VCF files
  + Slides: []()
  + Video: []()
+ Coffee Break and catch up (3:15-3:30)
+ Interpreting Exomiser output in excel/html (Video 10, 15 min: 3:30-4:00)
  + Download exomiser output in TSV format, open and explore tabs in excel
  + Visualize output HTML file
  + Slides: []()
  + Video: []()
+ *Problem Set 5* (4:00-5:00)

### Session 6: Wednesday, November 18th, 9:00AM-12:00PM (PST), SHY D308
+ *Review Problem Set 5* (9:00-10:00)
+ Open help session, working on Practical (10:00-12:00)

### Session 7: Wednesday November 18th, 3hr) CSB V2-222
+ *Practical* (2:00-5:00)
  + Each student presents their case to the group 
  + 15-20 minutes / presentation (10-15 min present, 5 min questions). 







