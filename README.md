# TrioCNV2
# Introduction 
TrioCNV2 is a tool designed to jointly detecting copy number variations from WGS data of the parent-offspring trio. TrioCNV2 first makes use of the read depth and discordant read pairs to infer approximate locations of copy number variations, and then employs the split read and local de novo assembly approach to refine the breakpoints.
# Installation
The easiest way to get TrioCNV2 is to download the binary distribution from the TrioCNV2 github release page. Alternatively, you can build TrioCNV2 from source with gradle.
1. git clone --recursive https://github.com/yongzhuang/TrioCNV2.git
2. Install gradle build tool (https://gradle.org/)
3. cd TrioCNV2 
4. gradle build 
You'll find the executable jar file in TrioCNV2/build/libs/. 

If you want to run TrioCNV2, you'll need:
1. Install Java SE Development Kit 8
2. Install R (Rscript exectuable must be on the path)
3. Install Runiversal (https://cran.r-project.org/web/packages/Runiversal/index.html) in R
4. Install String Graph Assembler (https://github.com/jts/sga) and add to the PATH

# Running
usage: java -jar TrioCNV2.jar [OPTIONS]  
**1. preprocess**  
   This option is used to extract the information from BAM files of the trio.

   usage: java -jar TrioCNV2.jar preprocess [OPTIONS]  

      -R,--reference    <FILE>   reference genome file (required)  
      -B,--bams         <FILE>   bam list file (required)  
      -P,--pedigree     <FILE>   pedigree file (required)  
      -M,--mappability  <FILE>   mappability file (required)  
      -O,--outputFile   <FILE>   output folder(required)  
         --deviation    <INT>    deletion insert size cutoff, median+deviation*SD(optional, default 6)  
         --window       <INT>    window size (optional, default 200)  
         --min_mapping_quality   <INT>    minumum mapping quality (optional,default 0)  
**2. call**  
   This option is used to jointly detect copy number variations.  

   usage: java -jar TrioCNV2.jar call [OPTIONS]  

      -I,--input           <FILE>   input folder got by the preprocess step (required)  
      -P,--pedigree        <FILE>   pedigree file (required)  
      -M,--mappability     <FILE>   mappability file (required)  
      -O,--output          <FILE>   output folder (required)  
         --exclude         <FILE>   exclude regions  
         --min_mappability <FLOAT>  minumum mappability(optional, default 0)  
         --mutation_rate   <FLOAT>  de novo mutation rate (optional, default 0.0001)  
         --transition_prob <FLOAT>  probability of transition between two different copy number states(optional, default 0.00001)  
         --min_distance    <INT>    minumum distance to merge two adjacent CNVs (optional, default 10K)  
         --outlier         <FLOAT>  the predefined percentage of outliers (optional, default 0.025)  
         --nt              <INT>    number of threads (optional, default 1)  
**3. refine**    
   This option is used to refine breakpoints.  

   usage: java -jar TrioCNV2.jar refine [OPTIONS]  

      -R,--reference    <FILE>   reference genome file (required)  
      -B,--bams         <FILE>   bam list file (required)  
      -P,--pedigree     <FILE>   pedigree file (required)  
      -I,--input        <FILE>   input folder got by the preprocess step (required)  
      -O,--output       <FILE>   output folder(required)  
         --deviation    <INT>    deletion insert size cutoff, median+deviation*SD(optional, default 3)  
         --size         <INT>    the size of expanded breakpoint regions (optional, default 400) 
         
# Input and Output

1. The 'preprocess' command will output a folder which contains the read depth,  discordant read pairs and estimated parameters.  

1) The bam list file is a single-column text file and each row is the path of one sample's BAM file in a trio.  
   Row 1: /path/to/father.bam  
   Row 2: /path/to/mother.bam  
   Row 3: /Path/to/offspring.bam  
   
2. The 'call' command will take the ouput folder of the first step and also output a folder which contains the candiate CNVs got by the read depth and discordant read pair approach.    

3. The 'refine' command will take the ouput folder of the second step and output a folder which contains the resulting SVs in the 'CNV.txt' file.   
   
   Column 1: chromosome  
   Column 2: start  
   Column 3: end  
   Coumnn 4: SAMPLE1's variant state   
   Column 5: SAMPLE2's variant state  
   Column 6: SAMPLE3's variant state
   Column 7: Evidence  
   

# Contact 
   yongzhuang.liu@hit.edu.cn
