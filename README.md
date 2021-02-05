# SimpactCyanPhylo

Extending Simpact Cyan to Phylogenetic Analysis.
This is a proof of the concept exercise, it shows how phylogenetic features capture HIV transmission dynamics. 
To reproduce the results, you need to follow instruction below and run the simulation. or load in your working environement the 'SimpactPaperPhyloExample' Rdata object which has all the results.


# SimpactCyanExamples

<!-- Created by David Niyukuri, 5 February 2020 -->


This readme content is part of what accompanies the article “SimpactCyan: An Open-source Simulator for Individual-Based Models in HIV Epidemiology with R and Python Interfaces”, available at: <https://github.com/wdelva/SimpactCyanExamples>


## CONTENTS

This repo contains the information necessary to reproduce the example:

* [Code and data files](#code-and-data-files)
   * Code file
   * Data files
* [System and software requirements](#system-and-software-requirements)
* [Copyright and licensing information](#copyright-and-licensing-information)
* [Contact information](#contact-information)

## CODE AND DATA FILES 

Three types of files (not including this readme file) are stored on this repository: code files, data files and figure files.


### Code file

All code is written in R. R is a statistical programming language and software package that is distributed under a GNU General Public License. R documentation and software is available for free download through the R Project for Statistical Computing website at http://www.r-project.org. The software is available as an executable file for a wide variety of operating systems and computer architectures, and a compilable binary is also available should you need it.

 
  SimpactCyanPhyloExample.R -- This file contains the script to generate the data, run the analysis and produce the figures for the HIV transmission network, time-resolved phylogenetic tree, and the distribution of internal nodes in the reconstructed phylogenetic tree and simulated HIV transmission events.


### Data files
  


  SimpactPaperPhyloExample.RData -- This file is an RData object, containing the data needed for the reproduction of the figures in the phylo example. These data are: transNet.yrs.Old (HIV transmission network object), dater.tree (time-stamped phylogenetic tree), i.vec (vector of calendar years), int.node.vec (vector of the number of internal nodes in one-year calendar time intervals), numC.tra (vector of the number of transmission events in one-year calendar time intervals), and pbtd (Viral lineage through time outputs from estimation of confidence intervals using parametric bootstrap for molecular clock dating).

  hiv.seq.C.pol.j.fasta -- This file is a fasta file of the HIV subtype C pol gene sequence of the seed infection (root of the phylogenetic tree).
  
  
### Figure files

  prevalence.pdf -- Age and gender stratified prevalence
  
  network.pdf -- Cumulative transmission network
  
  phylotree.pdf -- Phylogenetic tree built from generated viral sequence data
  
  events_nodes.pdf -- Transmission events and internal nodes of the phylogenetic tree
  
  LTT.pdf -- Plot for Viral Lineage Through Time
  

## SYSTEM AND SOFTWARE REQUIREMENTS

### Operating system

  We have only tested this code on personal computers (OS X Version 10.11.6 and Linux Ubuntu Version 16.04), and on the lengau cluster at the Cape Town Centre for High Performance Computing (CHPC) and the golett cluster of the Flemish Supercomputer Centre (VSC).

### Required software

  R version 3.4.4
  
  In order for the phylo example to work, you need to add the executable tools "Seq-Gen", and "FastTree" to your working directory, as well as the root viral gene sequence (hiv.seq.C.pol.j.fasta).

  Seq-Gen version 1.3.4. <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4> Simulates viral evolution across a transmission network and produces a sequence alignment. To install Seq-Gen, do the following:
  
  1. Visit the following Github repository to download the latest version of Seq-Gen: <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4>
  2. Click on the "Source Code" zip file to download
  3. Click on the zip file to unzip the folder
  4. Navigate to the source folder to confirm there is a file called "Makefile"
  5. Now you will need to compile the program using the Terminal on your computer
  6. Via the Terminal, change your working directory to the source folder by typing after the prompt: `cd "file/path/here/Seq-Gen-1.3.4 2/source"`
  7. Once your working directory has been set to the source folder, type after the prompt: `make`
  8. Now open the source folder and verify that a new file is present called "seq-gen"
  9. Copy that file and paste it into your R working directory

  FastTree version 2.1.10. <http://www.microbesonline.org/fasttree/#Install> Reconstructs a phylogenetic tree from a large alignment dataset. To install FastTree, do the following:
  
  1. Visit the website for downloading instructions: <http://www.microbesonline.org/fasttree/#Install>
  2. If you have a Linux operating system, you can directly download the executable files that are linked on that website. Those downloaded files can then be placed in your R working directory
  3. If you are using an OS X operating system, open the link "FastTree.c" in a new browser window
  4. Right-click on the program and click "Save as"
  5. Save anywhere on your computer
  6. Open the Terminal on your computer and change your working directory to the folder that contains "FastTree.c". After the prompt type:  `cd "file/path/here"`
  7. After the directory has been changed, after the prompt type: `gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm`
  8. Now check to see if a new executable file has been created in that folder
  9. Copy that file and paste it into your R working directory

  SimpactCyan version 0.21 and RSimpactCyan. SimpactCyan is the core program that allows fast simulation of HIV transmission across a sexual network. RSimpactCyan is the R package that enables initiation and running of models built by SimpactCyan. Installation instructions for both are at: <https://github.com/j0r1/RSimpactCyan/blob/master/INSTALLATION.md>

  A long list of auxiliary R packages is required to run the post-simulation analysis for the phylo example.

  install.packages("devtools")
  
  install.packages("pacman")
  
  library(devtools)

  install_github("j0r1/readcsvcolumns/pkg")

  install_github("wdelva/RSimpactHelp", dependencies = TRUE)

  p_load(Rcpp, ape, expoTree, data.table, readr, phangorn, dplyr, adephylo, treedater, geiger, picante, igraph, network, intergraph, ggtree, lubridate, ggplot2, ggnetwork, metafolio, magrittr, dplyr, tidyr, geomnet, RSimpactCyan, RSimpactHelper)
 

## COPYRIGHT AND LICENSING INFORMATION

All files are copyright protected and are made available under the GPL 3.0 License <https://www.gnu.org/licenses/gpl-3.0.en.html>. This means that this work is suitable for commercial use, that licensees can modify the work, that they must release the source alongside with Derivative Work, and that Derivative Work must be released under the same terms.


## CONTACT INFORMATION

Please contact Prof Wim Delva with questions regarding the files provided and their authorized usage:

David Niyukuri

Email: <NIYUKURI@sun.ac.za>



