########################################################
# Extracting HPDs etc. to a table, from a Beast MCC tree
########################################################

# Some of the functions in "tree_utils_v1.R", e.g. read_beast_prt, used for extracting the 
# bracketed statistics from BEAST NEXUS tree files (MCC files) are copied/lightly modified/heavily 
# modified from the R package " phyloch", by Christoph Heibl, licensed under GPL (>=2), available 
# at: http://www.christophheibl.de/Rpackages.html. Please cite phyloch also, also if you use this 
# feature.

# Load the R libraries for dealing with phylogenies & XML
library(XML)
library(ape)   # for read.tree
library(gdata) # for trim
library(XLConnect)    # for readWorksheetFromFile (seems to work better than read.xls)
library(devtools, httr) # for GitHub & Gist
library(BioGeoBEARS) # for sourceall

# Online source for BEASTmasteR (yes, I will make an official R package eventually, 
# but I am busy)
# Online (you can also download each and source locally)
library(devtools,httr)
source_gist(id="https://gist.github.com/nmatzke/b2dbf78532eca734881a", filename="sourceall_git.R")
sources = sourceall_git(repo="nmatzke/BEASTmasteR")

#########################################################
# BACKGROUND: the function 'read_beast_prt'
#########################################################
# read_beast_prt uses:
#
# - the BioGeoBEARS 'prt' (print tree to table) function, 
#   ...and combines it with...
# - 'extractBEASTstats2', heavily modified from phyloch
# 
# ...to produce a table (an R data.frame).
# 
# This table contains the tree in table format, 
# with each branch as a row.  The rows are in the same 
# order as the default node numbers in R -- that is, in the table:
# 
# - Rows 1-ntips are the tip nodes (and the branches below each tip node)
# - Row (ntips+1) is the root node. 
#      (Note: Typically, the root node is the only one that 
#             doesn't have a branch below it. It can have a branch, 
#             if root.edge is specified in the APE tree structure. 
#             Typically, MCC trees don't have this.)
# - The rest of the rows are internal nodes.
#

#########################################################
# SEEING THE 'APE' NODE NUMBERING SCHEME IN A PLOT
#########################################################
# See the first part of this BioGeoBEARS script:
#
# http://phylo.wikidot.com/example-biogeobears-scripts#node_numbering
# 

#########################################################
# EXAMPLE TREE FILE: 
#########################################################
#
# 'treeLog.mcc', in
# http://phylo.wikidot.com/local--files/beastmaster/LW12.zip
#
# ...or...
#
# http://phylo.wdfiles.com/local--files/beastmaster/dino_treeLog.mcc

# NEXUS tree file name (fn)
# E.g., an MCC tree from Beast/TreeAnnotator

# Example of a locally-stored MCC file
# nexfn = "treeLog.mcc"

# Note: check your working directory with
# getwd()
# list.files()

# Here, we'll use the remotely-stored MCC file
nexfn = "http://phylo.wdfiles.com/local--files/beastmaster/dino_treeLog.mcc"

# Read the tree to an APE tree object
tr = read.nexus(nexfn)

# Here's a tree table
trtable = prt(tr, printflag=FALSE)
head(trtable)

# What fields are there in the table?
names(prt)

# Here's the part of the table with the last tip nodes, the root node, and the
# beginning of internal nodes
trtable[87:93,]

# The tree has 89 tips and 88 internal nodes, so the table should have 89+88=177 nodes (rows)
dim(trtable)
# yep

###############################################
# EXTRACTING STATISTICS FROM A BEAST2 MCC TREE
###############################################
# Getting the tree table, with 
output = read_beast_prt(file=nexfn, digits = NULL, get_tipnames=TRUE, printflag=FALSE, include_rates=FALSE) 

# What's in output?
names(output)

# Here's the tree again
output$tr

# Here's the tree table, with stats added
head(output$prt_beast_nodestats)

# The tree has 89 tips and 88 internal nodes, so the stats table should have 89+88=177 nodes (rows)
prt_beast_nodestats = output$prt_beast_nodestats

# What fields are there?
names(prt_beast_nodestats)

# Q: But I don't see HPDs, or posterior probabilities,
#    when I use the head() command!
head(prt_beast_nodestats)

# A: TreeAnnotator only records HPD and posterior probabilities 
#    for internal branches/nodes.
# 
#    The posterior probability of a terminal branch (one immediately
#    below a tip node) is always 1, obviously.
#
#    The date of a tip node can vary (*if* you are doing tip-dating), 
#    but you can only access this tip-date and build a 95% HPD if you 
#    store the tip-dates in your Beast2 logfile, and then extract
#    them with some other BEASTmasteR functions. This is done 
#    in the example at: 
#    http://phylo.wikidot.com/beastmaster#tipdate_plots_script
# 

# See the stats appear for the internal nodes:
prt_beast_nodestats[87:93,]

###############################################
# BUT I WANT TO PLOT THE HPDs ON A TREE!
###############################################
# See this script / example:
# http://phylo.wikidot.com/beastmaster#tipdate_plots_script

###############################################
# Using the extracted data: example
###############################################

# Example plot: It is often the case that dating uncertainty
# increases as the median node age gets older. You can 
# check this with a command such as:

# internal node numbers
internal_nodenums = (length(tr$tip.label)+1) : (length(tr$tip.label)+tr$Nnode)
internal_nodenums

# Plot:
# x: time_bp (node time before present) 
# y: height_HPD_width (upper 95% cred. interval, minus lower 95% cred. interval)
#
# We are adding 66 my to the node ages, since the "top" of this tree is at the K-T boundary
plot(x=66+prt_beast_nodestats$time_bp[internal_nodenums], y=prt_beast_nodestats$height_HPD_width[internal_nodenums], xlab="node age (Ma)", ylab="HPD width", ylim=c(0, 50), xlim=c(250,65))
title("Theropod tree: Median node age vs. 95% HPD width")

# Uncertainty increases somewhat as we go back in time, but here the 
# effect is less than usual, probably because we have many fossils 
# throughout the study group

##############################################################################################
# IF YOU USE BEASTmasteR IN YOUR PUBLICATIONS, 
# PLEASE CITE IT: http://phylo.wikidot.com/beastmaster#citation
##############################################################################################