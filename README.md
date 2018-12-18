# GAGrank

Welcome to GAGrank! This software was developed in Joseph Zaia's lab in the Center for Biomedical Mass Spectrometry at the Boston University School of Medicine. GAGrank is a tool for determining glycosaminoglycan (GAG) structures using found tandem mass spectral fragments and bipartite network analysis. Currently, GAGrank is only available in command line form, but a web application is under development. Continue reading this file for instructions on using GAGfinder.

# Motivation

At the core of any sequencing method using MS<sup>2</sup> data is the relationship between the unknown sequence and its fragments: the actual sequence is ascertained based on the fragment ions generated in the fragmentation process. For GAGs, there are often many possible sequences for a given composition, and in an electron activated dissociation (ExD) experiment, there is a rich complement of product ions in the spectrum. The relationship between possible sequences and observed product ions is many-to-many, and can be represented in a network structure. In particular, the network structure is that of a bipartite network, which is a network whose nodes can be separated into two distinct partitions with edges only connecting nodes in one partition to nodes with the other partition. The concept of centrality in network analysis aims to solve this problem, and there are numerous existing algorithms for computing centrality measures. One such method is [PageRank](http://ilpubs.stanford.edu:8090/422/1/1999-66.pdf), developed by Brin and Page in 1996 for Google as a way to rank webpages according to their importance for search engine optimization purposes. Briefly, PageRank gives webpages higher importance values if they are linked to by other important webpages. PageRank was developed for general networks (i.e. not bipartite networks), but a recent method, [BiRank](https://arxiv.org/abs/1708.04396), was developed that adapts the PageRank algorithm for the specific case of bipartite networks. Briefly, BiRank gives nodes in partition A higher importance if they are linked to important nodes in partition B, and vice versa. Because of its design for bipartite networks, we employed BiRank with the goal of determining precursor sequence based on fragmentation patterns in the first GAG sequencing method developed using a network structure and network analysis algorithm, GAGrank. GAGrank was developed as a command line interface in the Python language.

# Installation

In progress.

# Running GAGrank command line interface

In progress.

# GAGrank arguments

Below are descriptions of all of the arguments for GAGrank:

```-h, --help  show this help message and exit
-c C        GAG class (required)
-i I        Input GAGfinder results file (required)
-r R        Reducing end derivatization (optional)
-m M        Precursor m/z (optional, but must be in mzML file)
-z Z        Precursor charge (optional, but must be in mzML file)
-s S        Number of sulfate losses to consider (optional, default 0)
-a A        Actual sequence, for testing purposes
```

# Example command

In progress.

# Output format

The output includes two self-explanatory columns: Sequence and GAGrank score. The sequences are string representations of the possible sequences for the given composition. The GAGrank score is the result of the BiRank algorithm on the network. GAGrank scores are not comparable across analyses, due to different numbers of possible sequences for different compositions.

# Conclusion

Any questions can be e-mailed to jzaia@bu.edu. Good luck!
