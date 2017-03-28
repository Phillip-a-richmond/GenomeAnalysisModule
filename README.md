# Introduction to Genomic Analysis Workshop Series
> This is the website for the UBC Advanced Research Computing Genomic Analysis Workshop Series.  Our mission is to teach genomic analysis and connect compute resources to researchers at UBC and across Canada. 

## Course Outline



Why? Well the 

## General Course Information & Prerequisites
The course is intended for academic researchers at Canadian institutions in Western Canada, that have access to nationally supported Research Computing through Compute Canada, specifically the WestGrid branch.  
+ **Prerequisites**
  + Install Vidyo, Instructions [Here](https://github.com/Phillip-a-richmond/ARC-Bioinformatics-Training/blob/master/Vidyo_instructions.md)
  + If you're a temporary user, refer to the email with instructions to log-in
  + For mac/linux users, only need native terminal which comes with the operating system.  You can access it in the Applications section for Mac users.  You also need to have X windows installed [Link](https://www.xquartz.org/)
  + For PC users, download and install [MobaXterm](http://mobaxterm.mobatek.net/) 
    + [Advanced MobaXterm usage](https://www.youtube.com/watch?v=Gkl8LD1rwlU) 
  + Initial experience in Linux Shell.  If you have NO experience in Linux, please visit another course of mine and go through modules 1-4 before coming to class. [Website](http://phillip-a-richmond.github.io/Bioinformatics-Introductory-Analysis-Course/)
  + IGV installed on local machine 
    + [IGV install](https://www.broadinstitute.org/software/igv/log-in) 
  + Filezilla installed on local machine
    + [FileZilla](https://filezilla-project.org/)
+ **Resources**
    + [Linux/Unix Cheatsheet](https://github.com/Phillip-a-richmond/ARC-Bioinformatics-Training/blob/master/UnixCheatSheet.pdf) 
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

2. Complete the Pre-class Quiz found [Here](https://github.com/Phillip-a-richmond/ARC-Bioinformatics-Training/blob/master/Pre-Course-Quiz.md)

3. Install Vidyo, IGV, FileZilla

## Course Details
### Session 1: Linux/Bioinformatics :
+ Monday, June 5th, 10:00-12:00
+ Linux - Logging in, filesystem, moving around, file editing, downloading files
+ Commands:
  + cp, ls, mv, cut, mkdir, rm, scp, wget
+ System tools: 
+ FTP/SCP client, terminal/shell
+ Bioinformatics  
+ File formats: BED, GFF, GTF, FASTA, DICT
+ Databases accessed: UCSC
  
### Session 2: Linux/Bioinformatics:
+ Tuesday, June 6th, 10:00-12:00
+ Linux - Shell scripts, interacting with the queue
+ Commands: 
  + sh, qsub, qstat, qdel, emacs, 
+ Bioinformatics  
+ File formats: FASTQ, FASTQ.GZ 
+ Tools used: FastQC 
+ Resources needed: ~1 core/person for 2 hour block  

### Session 3 (Short read mapping and visualization) (This is the piece I already taught)
+ Wednesday, June 7th, 10:00-12:00
+ Bioinformatics  
+ File formats: SAM, BAM, CRAM, BAI, sorted.BAM,
+ Tools used: BWA (with different options), samtools, IGV (downloaded on their own computer)
+ Datasets: DNA-sequencing (Human)
+ Resources needed: ~4 cores/person for 2 hour block  

### Session 4: Variant Calling - Small variants
+ Thursday, June 8th, 10:00-12:00
+ *Bioinformatics*  
+ File formats: VCF
+ Tools used: Picard, GATK, Platypus, Freebayes, Samtools mpileup
+ Datasets: DNA-sequencing (Human)
+ Resources needed: ~4 cores/person for 2 hour block  


### Mid-series quiz
+ Friday, June 9th


### Session 5: Variant Annotation
+ Monday, June 12th, 10:00-12:00    
+ Bioinformatics
+ File formats: VCF, BEDPE
+ Tools used: GEMINI
+ Datasets: DNA-sequencing (Human)
+ Resources needed: ~4 cores/person for 2 hour block  

### Session 6: RNA-seq
+ Tuesday, June 13th, 10:00-12:00
+ Bioinformatics
+ File formats: GTF review, 
+ Tools used: Tophat, cufflinks, cuffdiff (or updated version of tuxedo suite)
+ Datasets: RNA-sequencing (Mouse?, Human?, Worm, Fly)
+ Resources needed: ~4 cores/person for 2 hour block  

### Session 7: RNA-seq II
+ Wednesday, June 14th, 10:00-12:00
+ Bioinformatics
+ HTSeq → DESeq  





