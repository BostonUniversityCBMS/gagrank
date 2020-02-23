# GAGrank

Welcome to GAGrank! This software was developed in Joseph Zaia's lab in the Center for Biomedical Mass Spectrometry at the Boston University School of Medicine. GAGrank is a tool for determining glycosaminoglycan (GAG) structures using found tandem mass spectral fragments and bipartite network analysis. Currently, GAGrank is only available in command line form, but a web application is under development. Continue reading this file for instructions on using GAGrank.

# Motivation

At the core of any sequencing method using MS<sup>2</sup> data is the relationship between the unknown sequence and its fragments: the actual sequence is ascertained based on the fragment ions generated in the fragmentation process. For GAGs, there are often many possible sequences for a given composition, and in an electron activated dissociation (ExD) experiment, there is a rich complement of product ions in the spectrum. The relationship between possible sequences and observed product ions is many-to-many, and can be represented in a network structure. In particular, the network structure is that of a bipartite network, which is a network whose nodes can be separated into two distinct partitions with edges only connecting nodes in one partition to nodes with the other partition. The concept of centrality in network analysis aims to solve this problem, and there are numerous existing algorithms for computing centrality measures. One such method is [PageRank](http://ilpubs.stanford.edu:8090/422/1/1999-66.pdf), developed by Brin and Page in 1996 for Google as a way to rank webpages according to their importance for search engine optimization purposes. Briefly, PageRank gives webpages higher importance values if they are linked to by other important webpages. PageRank was developed for general networks (i.e. not bipartite networks), but a recent method, [BiRank](https://arxiv.org/abs/1708.04396), was developed that adapts the PageRank algorithm for the specific case of bipartite networks. Briefly, BiRank gives nodes in partition A higher importance if they are linked to important nodes in partition B, and vice versa. Because of its design for bipartite networks, we employed BiRank with the goal of determining precursor sequence based on fragmentation patterns in the first GAG sequencing method developed using a network structure and network analysis algorithm, GAGrank. GAGrank was developed as a command line interface in the Python language.

# Installation

Download your system specific version of GAGrank in the release tab. You will find an executable file in the `bin/dist` directory, which is GAGrank.

# Running GAGrank command line interface

In order to run GAGrank, you need to open a Command Prompt and make your way to the directory `gagrank\bin\dist`. If you are new to the Command Prompt, you can read a tutorial about different commands at (http://www.cs.princeton.edu/courses/archive/spr05/cos126/cmd-prompt.html). For instance, if you unzipped the `gagrank.zip` file in your `Downloads` folder, you could execute the following commands to get into the proper directory and list the files and folders (everything after the '>' symbol):

```[prompt]>C:
[prompt]>cd Users\[your username]\Downloads\gagrank\bin\dist
[prompt]>dir
```

There should be an executable file called `gagrank.exe`. This is your executable file to run GAGrank. Execute the following command to ensure that the installation was successful.

`[prompt]>gagrank.exe --help`

You should now see in your Command Prompt a list of arguments. These are important inputs from you that tell GAGrank exactly how to run.

# GAGrank arguments

Below are descriptions of all of the arguments for GAGrank:

```-h, --help  show this help message and exit
-c C        GAG class (required)
-i I        Input GAGfinder results file (required)
-r R        Reducing end derivatization (optional)
-m M        Precursor m/z (optional, but must be in mzML file)
-z Z        Precursor charge (optional, but must be in mzML file)
-s S        Number of sulfate losses to consider (optional, default 0)
-a A        Actual sequence, for testing purposes (optional)
```

# Example command

Now that you know where the executable file is and know about all of the arguments, let's go through an example analysis. Let's say you are analyzing an HS spectrum with precursor m/z of 362.4143 and precursor charge of -5. You ran the mass spectrometer in NETD mode. You already ran GAGfinder to determine the fragment ions in this spectrum, and the GAGfinder output file is located at `C:\Example\Not\Real\test_centroid.txt`. The reducing end tag adds an ethyl group (-CH2CH3). and you believe that there may be up to two sulfate losses. You would use the following command (everything after the '>' symbol):

`[prompt]>gagrank.exe -c HS -i C:\Example\Not\Real\test_centroid.tsv -r CH2CH3 -m 362.4143 -z -5 -s 2`

The output would then be found at `C:\Example\Not\Real\test_GAGrank_results.tsv`.

# Output format

The output includes two self-explanatory columns: Sequence and GAGrank score. The sequences are string representations of the possible sequences for the given composition. The GAGrank score is the result of the BiRank algorithm on the network. GAGrank scores are not comparable across analyses, due to different numbers of possible sequences for different compositions.

# Conclusion

Any questions can be e-mailed to jzaia@bu.edu. Good luck!
