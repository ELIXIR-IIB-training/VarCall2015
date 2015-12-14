---
course: NGS for evolutionary biologists: from basic scripting to variant calling
title: Tutorial - Variant calling, filtering - exercises
requires:
author:  Chiara Batini, Pille Hallast
credits: Matthew Blades, Kate Lee (BBASH, University of Leicester, UK)
time:
---
------------
> #### Learning Objectives
---


# Getting Started with this handbook


Each stage of the workflow will be described in terms of what is being done and the command required to run the software.  Commands will be in the following format:

```
>$ command you need to type to run the software
```

Where:


`>$` is the command prompt. The actual command prompt on your terminal will be in a format something like this:

```
your_user_name@node013.pico:[~]$
```

e.g in my case the command prompt looks like this:
```
phallast@node013.pico:[~]$
```

However, to save space in this manual I have replaced the actual command prompt with this:
`>$` or removed at all.

The text in gray areas is the actual text you need to type into the command window to run the software.
**NOTE:** After typing in a command, hit return to execute that command


A command is generally the name of the software to be used followed by a list of parameters and options for how the software should be run, e.g. the command below will run the software *cutadapt*. The text after cutadapt i.e.  `-a`, `-i` and `-o` are all *cutadapt* parameters followed by their input values:

```
>$ cutadapt -a ATGAATCTA -i data_in.fastq -o data_out.fastq
```

**NOTE:** Be aware of spaces between text in commands as they are important!


During the practicals questions will be asked and a blank space left for you to fill in your answers.   Questions are written in green text.
