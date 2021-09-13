# start with the uncompressed sequence reads (.fq)

# run HISAT2 to align seqence reads to reference genome
hisat2
    -q  # (query files are in .fq format)
    -p 16  # (use 16 threads for alignment)
    --pen-noncansplice 1000000  # (penalty for a non-canonical splice site)
    -x path/to/index  # (Index filename prefix (minus trailing .X.ht2))
    -1 input.fq  # (Files with #1 mates, paired with files in -2)
    -2 input.fq  # (Files with #2 mates, paired with files in -1)
    -S output.sam  # (File for SAM output)

# log the alignment results, importantly the alignment rate (should be >70%)

# after this point, .fq files can be gzipped again to save space

# convert .sam files to more compact binary .bam files
samtools view
    -b --threads 16  # (output in binary .bam format, use 16 threads for compression)
    input.sam > output.bam  # (input .sam file and output .bam file)

# after this point, .sam files can be deleted to save space

# sort the .bam files using samtools sort
samtools sort
    --threads 16  # (use 16 threads for sorting)
    input.bam -o output.sort.bam  # (input .bam and output sorted .bam)

# after this step we can delete the unsorted .bam files to save space

# count reads to genomic features
featureCounts
    -p  # (fragments or templates will be counted instead of reads) 
    -t exon  # (specify feature type in GTF annotation)
    -a gencode.VM25.annotation.gtf  # (name of annotation file, GTF format)
    -g gene_name  # (specify attribute type in GTF annotation)
    -T 16  # (number of threads)
    -o output.txt  # (name of output file including read counts)
    f1.sort.bam f2.sort.bam ... fN.sort.bam  # (list all sorted .bam files to use)

