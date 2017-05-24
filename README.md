# Introduction to Genomic Analysis Workshop Series
> This course is intended for Canadian researchers interested in learning the basic steps of genomic analysis using high performance computing. Using the computing resources of Compute Canada, specifically the WestGrid region, participants will work through practical exercises in the workshop sessions and receive a certificate of completion at the end of the series.  



---

## General Course Information & Prerequisites
### Registering for the course
Please register for the course using our EventBrite Link [here](https://www.eventbrite.ca/e/introduction-to-genomic-analysis-workshop-series-registration-33475860199)

### Course Format

> The course is intended for academic researchers at Canadian institutions in Western Canada, that have access to nationally supported Research Computing through Compute Canada, specifically the WestGrid branch.  

+ 2 hour workshops:
  + 45 minute lecture
  + 15 minute Q&A
  + 45 minute problemset
  + 15 wrap up
+ Each workshop has a prerequisite "quiz".  The problemset from the preceding workshop is the "quiz"
+ If you're familiar with a given subject, you just need to be able to complete the entry quiz instead of attending that subject
  + *Example 1*: Familiar with Linux, but want to learn about short read mapping, do the quiz for Session 2, and show up for Session 3. 
  + *Example 2*: Familiar with short read sequencing analysis and visualization, but want to learn about variant calling, do the quiz for Session 3 and show up to Session 4.
+ Each workshop will have local environment, hosted at BC Children's Hospital with refreshments provided
+ Each workshop will also have a virtual environment, hosted through Vidyo (See details on Vidyo [below](#Software-Installs) )

### Prerequisites  

#### Getting a Compute Canada WestGrid Account  
+ Follow the instructions [here](https://www.computecanada.ca/research-portal/account-management/apply-for-an-account/)
+ If you're a temporary user, refer to the email with instructions to log-in

#### Getting a temporary Compute Canada WestGrid Account  
+ IF YOU ARE UNABLE TO DO THE ABOVE!!!
+ Then a temporary account will be provided to you in the week leading up to the course, and will dissolve when the course completes.

#### Software Installs  
+ Install Vidyo, Instructions [Here](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Resources/Vidyo_instructions.md)
+ For mac/linux users, only need native terminal which comes with the operating system.  You can access it in the Applications section for Mac users.  You also need to have X windows installed [Link](https://www.xquartz.org/)
+ For PC users, download and install [MobaXterm](http://mobaxterm.mobatek.net/) 
  + [Advanced MobaXterm usage](https://www.youtube.com/watch?v=Gkl8LD1rwlU)  
+ IGV installed on local machine 
  + [IGV install](https://www.broadinstitute.org/software/igv/log-in) 
+ Filezilla (or other file transfer client) installed on local machine
  + [FileZilla](https://filezilla-project.org/)
  
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
+ WestGrid resources 
  + [WestGrid website](https://www.westgrid.ca/) 
  + [Running jobs](https://www.westgrid.ca/support/running_jobs) 
  + [Software available](https://www.westgrid.ca/support/software/) 
      
### To-Do before the first class  
1. Get your login information for WestGrid, and log-in to the orcinus compute server.
ssh \<username\>@orcinus.westgrid.ca
example:
ssh richmonp@orcinus.westgrid.ca

2. Complete the Pre-class Quiz found [Here](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/PreQuiz_Session1.txt)

3. Install Vidyo, IGV, FileZilla

---

## Course Details  

### Session 1: Introduction to Linux I, Basic Command Line
> Monday, June 5th, 10:00AM-12:00PM (PST)
+ [Session 1 Pre-Quiz](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/PreQuiz_Session1.txt)
+ [Session 1 Slides](https://docs.google.com/presentation/d/1vaHO9tewJhnpn3CQkIADKnJka4SrKCFb2_zOMTOi7yc/edit#slide=id.g1957375280_1_15)
+ *Linux* 
  + Logging in, filesystem, moving around, file editing, downloading files
  + Commands:
  + cp, ls, mv, cut, mkdir, rm, scp, wget
+ *System tools*: 
  + terminal/shell
+ *Bioinformatics*  
  + File formats: BED, GFF/GTF, VCF
+ [Github Gist](https://gist.github.com/Phillip-a-richmond/8b36d58fb2f03d4760869b0676d50692)
+ [Session 1 ProblemSet](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session1_Problemset.txt)

---

### Session 2: Introduction to Linux II, Interacting with the Queue
> Tuesday, June 6th, 10:00AM-12:00PM (PST)
+ [Session 2 Pre-Quiz](https://github.com/Phillip-a-richmond/ARC-Bioinformatics-Training/blob/master/Pre-Course-Quiz.md)
+ [Session 2 Slides](https://docs.google.com/presentation/d/1u5R23969caKTCPB3YIYRZRNW2zIcktZJC5S89D-GgjU/edit#slide=id.g1deaa80ad7_0_0)
+ *Linux* 
  + File editor
  + FTP client
  + Shell scripts, interacting with the queue
  + Commands: sh, qsub, qstat, qdel, emacs, 
+ *Bioinformatics* 
  + File formats: FASTQ, FASTQ.GZ 
  + Tools used: FastQC 
  + HPC Resources needed: ~1 core/person for 2 hour block  
+ [Github Gist](https://gist.github.com/Phillip-a-richmond/76baff0d4525817e0ca59a638bb374a0)
+ [Session 2 ProblemSet](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session2_Problemset.txt)

---

### Session 3: Short read mapping and visualization
> Wednesday, June 7th, 9:30AM-11:30AM (PST)
+ [Session 3 Pre-Quiz](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/Assignments/Session2_Problemset.txt)
+ [Session 3 Slides](https://docs.google.com/presentation/d/1u5R23969caKTCPB3YIYRZRNW2zIcktZJC5S89D-GgjU/edit#slide=id.g1deaa80ad7_0_0)
+ *Bioinformatics*  
  + File formats: SAM, BAM, CRAM, BAI, sorted.BAM,
  + Tools used: BWA (with different options), samtools, IGV (downloaded on own computer)
  + Datasets: DNA-sequencing (Human)
  + HPC Resources needed: ~4 cores/person for 2 hour block
+ [Session 3 ProblemSet]()

---

### Session 4: Variant Calling - Small variants
> Thursday, June 8th, 10:00AM-12:00PM (PST) 
+ [Session 4 Pre-Quiz](https://github.com/Phillip-a-richmond/ARC-Bioinformatics-Training/blob/master/Pre-Course-Quiz.md)
+ [Session 4 Slides]()
+ *Bioinformatics*  
  + File formats: VCF
  + Tools used: Picard, GATK, Platypus, Freebayes, Samtools mpileup
  + Datasets: DNA-sequencing (Human)
  + HPC Resources needed: ~4 cores/person for 2 hour block  
+ [Session 4 ProblemSet]()

---

### Mid-series Exam
> Friday, June 9th
+ [Mid-Series Exam]()

---

### Session 5: Variant Annotation
> Monday, June 12th, 10:00AM-12:00PM (PST)   
+ [Session 5 Pre-Quiz](https://github.com/Phillip-a-richmond/ARC-Bioinformatics-Training/blob/master/Pre-Course-Quiz.md)
+ [Session 5 Slides]()
+ *Bioinformatics*
  + File formats: VCF, BEDPE
  + Tools used: GEMINI
  + Datasets: DNA-sequencing (Human)
  + HPC Resources needed: ~4 cores/person for 2 hour block  
+ [Session 5 ProblemSet]()

---

### Session 6: RNA-seq I, Read Mapping
> Tuesday, June 13th, 10:00AM-12:00PM (PST) 
+ [Session 6 Pre-Quiz](https://github.com/Phillip-a-richmond/ARC-Bioinformatics-Training/blob/master/Pre-Course-Quiz.md)
+ [Session 6 Slides]()
+ *Bioinformatics*
  + File formats: GTF review
  + Tools used: Tophat, cufflinks, cuffdiff (or updated version of tuxedo suite)
  + Datasets: RNA-sequencing (Mouse?, Human?, Worm, Fly)
  + HPC Resources needed: ~4 cores/person for 2 hour block  
+ [Session 6 ProblemSet]()

---

### Session 7: RNA-seq II, Differential Expression
> Wednesday, June 14th, 9:30AM-11:30AM (PST)
+ [Session 7 Pre-Quiz](https://github.com/Phillip-a-richmond/ARC-Bioinformatics-Training/blob/master/Pre-Course-Quiz.md)
+ [Session 6 Slides]()
+ *Bioinformatics*
  + HTSeq → DESeq  
+ [Session 7 ProblemSet]()



## Proudly Supported By:
+ University of British Columbia Advanced Research Computing
+ BC Children's Hospital Evidence 2 Innovation
+ Compute Canada WestGrid

![BCCH Logo](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/bcch_logo1.png)

![WestGrid Logo](https://github.com/Phillip-a-richmond/Introduction-to-Genomic-Analysis/blob/master/wesgrid_logo_2016.png)





