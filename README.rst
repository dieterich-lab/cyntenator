==============================================
Cyntenator : Progressive Gene Order Alignments
==============================================

Cyntenator is a software for identification of conserved syntenic blocks between multiple genomes. 
The program computes Smith-Waterman alignments of sequences, whereby the alphabet consists of all 
annotated genes and the scoring system is defined by protein sequence similarities and distances 
between species in a phylogenetic tree. The algorithm is an extension of the Syntenator partial order 
aligner, described in Rödelsperger and Dieterich, 2008. 

Content:
""""""""""""""""""

1. Installation
2. Usage
3. File Format
4. License


***************
1. Installation
***************

In order to install cyntenator on linux, simply type

.. code:: bash

 $ make


this will create a binary executable named cyntenator. The Makefile uses the g++ compiler, 
in general it should be possible to compile the program on windows but this was not tested.

********
2. Usage
********

If cyntenator is called without any arguments, a help message is printed on standard output

.. code:: bash

    ./cyntenator

    program -t guide-tree -h homology_type ... 
    
    guide-tree:
            -t "((rat.txt mouse.txt ) human.txt)"
    Homology:
            -h id
            -h blast [file]
            -h orthologs [file]
            -h phylo [file] [weighted_tree]
    
    Alignment Parameters:
            -thr    threshold (4)
            -gap    gap (-2)
            -mis   mismatch (-3)
    
    Filter options:
            -filter [int] best alignments or only unique assignments n=0 (100)
            -coverage [int] each gene may occur only c times in alignments (2)
            -length [int] minimum alignment length treshold (1)
            -last prints only the alignments at the last step
    
    Output:
            -o output file



Obligatory arguments are the guiding tree which specifies the order in which genomes are aligned. It can be 
specified by a notation in parentheses (e.g. "((rat.txt mouse.txt ) human.txt)"). In this example rat.txt, mouse.txt, and human.txt are the gene annotation files for the respective species.

The second argument which is needed is the type of homology and a corresponding similarity file. The type
"id" assigns a match, whenever two IDs (GENE Symbols, Transcription factor Binding sites or any other Marker)
are equal. The "orthology" type assigns matches for pair of genes that are given in a correspondence file, furthermore
correspondences can be assigned similarity scores like BLAST bitscores, this is the option, that we used in 
the Syntenator approach.  In order to weight matches by sequence similarity and phylogenetic distance,
in addition to the BLAST file, a weighted species tree has to be defined with weights specified by [:WEIGHT] 
after each closing bracket (e.g. "((HSX.txt:1.2 MMX.txt:1.3):0.5 CFX.txt:2.5):1" ). The final weight at the root 
note can be omitted.

Alignment parameters can be modified via the options thr gap and miss. 


***************
3. File Formats
***************

There are three possible input formats which can be used as annotation files. The format has to be specified in 
the first row of the file:

#genome
#sequence
#alignment

#genome this is the standard format for alignment of genomes. Annotations can be obtained from the ensembl database
and should contain white space separated gene coordinates (ID,CHROMOSOME,START,END,STRAND, see HSX.txt as example). 
The same is true for the #sequence format which should be used for alignment using IDs as similarity measure. 
The only difference is that the genome format requires unique IDs and merges redundant gene annotations which come 
from different transcripts.

The alignment format is the output format of cyntenator. In the header the original annotation files
are enlisted. The following lines show all computed alignments, sorted by score and uniqueness.
An Alignment starts with the word Alignment followed by a number and then a score. Each following line
starts with a score for this row and the aligned genes with their orientation. "-" characters denote gaps.

.. code:: bash

    #alignment CFX.txt HSX.txt MMX.txt

    Alignment 1 9.18021
    0.348995  ENSCAFG00000016080  (+)  ENSG00000184205  (+)  ENSMUSG00000041096  (+)  
    1.89445   ENSCAFG00000016133  (-)  ENSG00000126012  (-)  ENSMUSG00000025332  (-)  
    0.348995  -                   (.)  ENSG00000204349  (+)  -                   (.)  
    1.36421   ENSCAFG00000016179  (-)  ENSG00000124313  (-)  ENSMUSG00000041115  (-)  
    2.88981   ENSCAFG00000016241  (-)  ENSG00000072501  (-)  ENSMUSG00000041133  (-)  
    3.3294    ENSCAFG00000016274  (+)  ENSG00000158423  (+)  ENSMUSG00000025257  (+)  
    3.77215   ENSCAFG00000016277  (-)  ENSG00000072506  (-)  ENSMUSG00000025260  (-)  
    5.3176    ENSCAFG00000016311  (-)  ENSG00000086758  (-)  ENSMUSG00000025261  (-)  
    6.71358   ENSCAFG00000016389  (-)  ENSG00000172943  (-)  ENSMUSG00000041229  (-)  
    5.16812   -                   (.)  -                (.)  ENSMUSG00000067230  (+)  
    6.57727   ENSCAFG00000016414  (-)  ENSG00000184083  (-)  ENSMUSG00000025262  (-)  
    5.03181   -                   (.)  -                (.)  ENSMUSG00000072901  (+)  
    6.57727   ENSCAFG00000016427  (-)  ENSG00000196632  (-)  ENSMUSG00000041245  (-)  
    6.86235   ENSCAFG00000016440  (+)  ENSG00000158526  (+)  ENSMUSG00000025264  (+)  
    8.19742   ENSCAFG00000016444  (-)  ENSG00000102302  (-)  ENSMUSG00000025265  (-)  
    9.18021   ENSCAFG00000016455  (+)  ENSG00000130119  (+)  ENSMUSG00000025266  (+)  


It is possible to compute pairwise alignments parallel and then compute the multiple genome alignments
by specifying the alignment output files as input.  The following 2 runs are test runs to see, whether 
cyntenator can read pairwise alignments

PART 1: compute pairwise, read pairwise, compute triple
=======================================================

.. code:: bash

 ./cyntenator -t "(HSX.txt MMX.txt)"     -h phylo HSCFMM.blast  "((HSX.txt:1.2 MMX.txt:1.3):0.5 CFX.txt:2.5):1" > human_mouse
 ./cyntenator -t "(human_mouse CFX.txt)" -h phylo HSCFMM.blast  "((HSX.txt:1.2 MMX.txt:1.3):0.5 CFX.txt:2.5):1" > test1.txt 


PART 2: compute pairwise, compute triple
=======================================================

.. code:: bash

 ./cyntenator -last  -t "((HSX.txt MMX.txt) CFX.txt)" -h phylo HSCFMM.blast  "((HSX.txt:1.2 MMX.txt:1.3):0.5 CFX.txt:2.5):1" > test2.txt

test1.txt and test2.txt will contain the same alignments in the end. Not, that the -last option is used in the second 
example in order to report only the alignments created from the HSX-MMX vs CFX multiple species comparison.


***************
4. License
***************

This software was written by Christian Rödelsperger  and is distributed under GNU General Public License (GPL)


