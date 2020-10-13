# Genome Analysis Module
> This course is intended for trainees in the Genetic Counselling MSc Program at the University of British Columbia. The computing resources are hosted at BC Children's Hospital Research Institute. All data is synthetic, and for more information contact Dr. Phillip Richmond (prichmond@bcchr.ca). If you are looking for content related to the manuscript, please see Release 1. 

---
## License
> This work is distributed under the creative commons license (CC BY-SA 3.0).  Details found [here](https://creativecommons.org/licenses/by-sa/3.0/deed.en_US)

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
+ *Problem Set 1* (4:00-5:00)
  + Link to Problem Set: []()

### Session 2: Wednesday, November 4th, 9:00AM-12:00PM (PST)
+ *Problem Set 1 Review* (9:00-10:00)
+ Interacting with the cluster scheduler (Lecture: 10:00-10:30)
  + Learning the cluster
  + SLURM job submitting
  + Module system for software
+ Interacting with the cluster scheduler (Video 4, 15 min: 10:30-11:00)
  + Edit job scripts, submit to queue, collect output
+ Coffee Break and catch up (11:00-11:15)
+ *Problem Set 2* (11:15-12:00)

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

### Session 4: Monday, November 9th, 1:00PM-5:00PM (PST)
+ Problem Set 3 Review (1:00-2:00)
+ Variant calling and visualization (Lecture: 2:00-2:30)
  + Many tools for variant calling, but require same inputs
  + Different classes of variants which can be called 
  + Call variants with DeepVariant (Video 6, 15 min: 2:30-3:00) 
  + Call variants using pre-built DeepVariant scripts on pre-processed BAM files
+ Coffee Break and catch up (3:00-3:15)
+ Visualization with IGV (Video 7, 10 min: 3:15-4:00)
  + Download and install IGV
  + Load GRCh38 genome
  + Load read files
  + Zoom to region with read coverage
  + Take snapshot of variant
+ Problem Set 4 (4:00-5:00)



_*Any assistance needed with Problem Sets 1-4 covered in the week between sessions. At this point, we need to ensure that all students have processed their samples to the VCF stage.*_










+ *Linux* 
  + Logging in, filesystem, moving around, file editing, downloading files
  + Commands:
  + cp, ls, mv, cut, mkdir, rm, scp, wget
+ *System tools*: 
  + terminal/shell
+ *Bioinformatics*  
  + File formats: BED, GFF/GTF, VCF
+ [Github Gist](https://gist.github.com/Phillip-a-richmond/8b36d58fb2f03d4760869b0676d50692)
+ [Session 1 Etherpad](https://etherpad.openstack.org/p/Introduction_To_Genomic_Analysis_Session1)
+ [Session 1 ProblemSet](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session1_Problemset.txt)
+ [Link to raw recording](https://vidyoreplay.computecanada.ca/replay/showRecordingExternal.html?key=bBtILur7bxNhjMp)
  + Note: Fast-forward to 6:33 for presentation to start
+ [Link to final recording]()


---

### Session 2: Introduction to Linux II, Interacting with the Queue
> Tuesday, June 6th, 10:00AM-12:00PM (PST)
+ [Session 2 Pre-Quiz](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session1_Problemset.txt)
+ [Session 2 Slides](https://docs.google.com/presentation/d/1u5R23969caKTCPB3YIYRZRNW2zIcktZJC5S89D-GgjU/edit#slide=id.g1deaa80ad7_0_0)
+ *Linux* 
  + File editor
  + FTP client
  + Shell scripts, interacting with the queue
  + Commands: sh, qsub, qstat, qdel, emacs, 
+ *Bioinformatics* 
  + File formats: FASTQ, FASTQ.GZ 
  + HPC Resources needed: ~1 core/person for 2 hour block  
+ [Github Gist](https://gist.github.com/Phillip-a-richmond/76baff0d4525817e0ca59a638bb374a0)
+ [Session2 Etherpad](https://etherpad.openstack.org/p/Introduction_To_Genomic_Analysis_Session2)
+ [Session 2 ProblemSet](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session2_Problemset.txt)
+ [Link to raw recording](https://vidyoreplay.computecanada.ca/replay/showRecordingExternal.html?key=invabc4OKUOtXxa)
  + Note: Fast-forward to 1:52 for presentation to start
+ [Link to final recording]()
  


---

### Session 3: Short read mapping and visualization
> Wednesday, June 7th, 9:30AM-11:30AM (PST)
+ [Session 3 Pre-Quiz](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session2_Problemset.txt)
+ [Session 3 Slides](https://docs.google.com/presentation/d/1BJJjY3TaR-TEe8HhBCu1pXp2HFCKgTraRkayrUdvG4Y/edit#slide=id.g1deaa80ad7_0_0)
+ *Bioinformatics*  
  + File formats: SAM, BAM, BAI, sorted.BAM,
  + Tools used: BWA (with different options), samtools, IGV (downloaded on own computer)
  + Datasets: DNA-sequencing (Human)
  + HPC Resources needed: ~4 cores/person for 2 hour block
+ [Github Gist](https://gist.github.com/Phillip-a-richmond/bf669287f33b5bd91c27549b94250834)
+ [Session 3 ProblemSet](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session3_Problemset.txt)
+ [Session 3 Etherpad](https://etherpad.openstack.org/p/Introduction_To_Genomic_Analysis_Session3)
+ [Link to raw recording](https://vidyoreplay.computecanada.ca/replay/showRecordingExternal.html?key=ycoB7nogbzccDOO)
  + Note: Fast-forward to 1:28 for presentation to start
+ [Link to final recording]()

---

### Session 4: Variant Calling - Small variants
> Thursday, June 8th, 9:30AM-11:30PM (PST) 
+ [Session 4 Pre-Quiz](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session3_Problemset.txt)
+ [Session 4 Slides](https://docs.google.com/presentation/d/1vdX5-nBWnxlQ2C09i8J1OgOwirj-EA-nLwvXOHR58ek/edit#slide=id.g1e4523bd13_1_30)
+ *Bioinformatics*  
  + File formats: VCF
  + Tools used: Picard, GATK, Platypus, Freebayes, Samtools mpileup
  + Datasets: DNA-sequencing (Human)
  + HPC Resources needed: ~4 cores/person for 2 hour block  
+ [Github Gist](https://gist.github.com/Phillip-a-richmond/4f6da08a09881cd088a49d49656ae243)
+ [Session 4 Etherpad](https://etherpad.openstack.org/p/Introduction_To_Genomic_Analysis_Session4)
+ [Session 4 ProblemSet](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session4_Problemset.txt)
+ [Link to raw recording](https://vidyoreplay.computecanada.ca/replay/showRecordingExternal.html?key=QJpWP4xxkhLwIyD)
  + Note: Fast-forward to 01:20 for presentation to start
+ [Link to final recording]()

---

### Mid-series Exam
> Friday, June 9th
+ Will be available for answering questions from 10:00-12:00
  + Virtually:  WG_Training room
  + Locally in BCCHR Rm# 2108
+ [Mid-Series Exam](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session4_Problemset.txt)
+ [Exam etherpad help](https://etherpad.openstack.org/p/Introduction_To_Genomic_Analysis_Session4)

---

### Session 5: Variant Annotation
> Monday, June 12th, 10:00AM-12:00PM (PST)   
+ [Session 5 Pre-Quiz](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/tree/master/Assignments)
+ [Session 5 Slides](https://docs.google.com/presentation/d/1WJ6K-uKoKKcpIA0joufkU4Ym5rF6w3wAahzPWWjHSMs/edit?usp=sharing)
+ *Bioinformatics*
  + File formats: VCF
  + Tools used: GEMINI
  + Datasets: DNA-sequencing (Human)
  + HPC Resources needed: ~4 cores/person for 2 hour block  
+ [Session 5 ProblemSet](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session5_Problemset.txt)
+ [Session 5 Etherpad](https://etherpad.openstack.org/p/Introduction_To_Genomic_Analysis_Session5)
+ [Session 5 Github Gist](https://gist.github.com/Phillip-a-richmond/7cce5f287ecc678d29acea980a3086b6)
+ [Link to raw recording](https://vidyoreplay.computecanada.ca/replay/showRecordingExternal.html?key=E2wWUWpLkAZE5aJ)
  + Note: Recording starts at 0:00.


---

### Session 6: RNA-seq I, Read Mapping, Transcript Quantification
> Tuesday, June 13th, 10:00AM-12:00PM (PST) 
+ [Session 6 Pre-Quiz](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session4_Problemset.txt)
+ [Session 6 Slides](https://docs.google.com/presentation/d/1QyHtpR6rSRqB_EPY7UuoQY_ufZ8Qp1pctgEQ2hflERs)
+ *Bioinformatics*
  + File formats: GTF review
  + Tools used: HISAT2, stringtie
  + Datasets: RNA-sequencing Inbred Mouse Brain Expression Data
  + HPC Resources needed: ~4 cores/person for 2 hour block 
+ [Session 6 Etherpad](https://etherpad.openstack.org/p/Introduction_To_Genomic_Analysis_Session5)
+ [Session 6 Github Gist](https://gist.github.com/Phillip-a-richmond/ec216db03b1e6c4b1f645934a5c6ebfe)
+ [Session 6 ProblemSet](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session6_Problemset.txt)
+ [Link to raw recording](https://vidyoreplay.computecanada.ca/replay/showRecordingExternal.html?key=Bb6vMbJqnF9iwqD)
  + Note: Skip ahead to 00:50 for presentation to start

---

### Session 7: RNA-seq II, Differential Expression
> Wednesday, June 14th, 9:30AM-11:30AM (PST)
