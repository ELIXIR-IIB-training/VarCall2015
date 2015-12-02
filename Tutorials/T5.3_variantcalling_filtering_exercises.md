---
course: NGS for evolutionary biologists: from basic scripting to variant calling
title: Tutorial - Variant calling, filtering - exercises
requires:
author:  Chiara Batini, Pille Hallast  
time:
---
------------
> #### Learning Objectives
------------


# Variant Calling and filtering practicals

**A starter:**

Is this an haploid or a diploid genome? Have you wondered?
How can you answer that question having a vcf file available?


**Now, imagine you have just got the final bam file from your first NGS run in your project. How would you deal with it?**

## Practical 1 – compare different SNPs variant callers

One way to go is to try different variant callers, and then keep the overlapping variants. As mentioned, GATK is another popular software for variant calling of SNPs and small indels. GATK implements two callers:
- UnifiedGenotyper (UG)
- HaplotypeCaller (HC)

Check the [GATK](https://www.broadinstitute.org/gatk/guide/tooldocs/) website to find out the difference among the two and decide which one is the best to use.

Compare the results of GATK with those of samtools using [vcftools](https://vcftools.github.io/perl_module.html#vcf-compare). Use the same final bam file you have used with samtools.

If you still have time, you might want to try out two additional callers that are available on pico:
- [FreeBayes](https://github.com/ekg/freebayes#readme)
- [VarScan2](http://dkoboldt.github.io/varscan/using-varscan.html)

## Practical 2 – apply a set of filters to your final vcf

Decided which caller to use, you might want to know how different filters affect your final variant list.

Start with this list of filters. Think about the order in which you should use them. Obviously this is not the right one ;)

- Filter indels out – better doing this during the variant calling (like we did with BQ and MQ) or after on the vcf file?
- DP
- Filter out monomorphic sites – why could they be there? They can…
- Check --max-mac option in vcftools to filter them out
- SnpGap
- Proximity to Indels
- Anything else you think is important among the filters available in vcf-annotate?


Have you wondered what would happen if you filtered for a different depth? You can do it and then use vcftools to [compare](https://vcftools.github.io/perl_module.html#vcf-compare) the two vcf files.


## Practical 3 – use Pindel to call structural variants (SVs) from your data

Now try to call SVs from your samples.

Small indels have also been called by samtools and GATK callers - try to compare the calls from these softwares. Can you find overlapping ones?

Let our aim be finding medium and large (>10bp) SVs. How would you filter the calls? How many insertions and and deletions are you left with?


If you still have time, there is an other popular SV calling software available on pico:
[BreakDancer](http://gmt.genome.wustl.edu/packages/breakdancer/index.html)

Try to understand what are the main differences between Pindel and BreakDancer. Try running BreakDancer on your own.
