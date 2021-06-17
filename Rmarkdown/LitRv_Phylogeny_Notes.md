## Phylogeny

* Three Different methods to generate trees 1) Genetic pairwise hamming distance 2) ML 3) Bayesian tree

* Distance tree are not reliable but fast to generate. Not usually accepted in journals.

* Maximum Likelihood phylogeny tree was generated using Mega software v7.0.26. The Tamura-Nei model with discrete Gamma distribution was applied to model evolutionary rate differences among sites (4 categories (+G, parameter = 0.1133))

* Whether or not we are generating ML or Bayesian tree we need a model to construct phylogeny. For ML tree, the best way to do this is to use `JModelTest`, in which you simply input your alignment and it will tell you the best model for your data (MEGA/PhyML/RAxML). But for Bayesian BEAST tree, I used `bmodeltest` which is a new method the site model can be inferred (and marginalized) during the MCMC analysis and **does not** need to be pre-determined, as is now often the case in practice, by likelihood-based methods (Mr.Bayes can also construct Bayesian tree).

* Tree Generation steps:
  * SNPs (alignment fasta)
  * Use a model (for example: GTR - an exact model depends on the data so we have to determine it using (j/b)modeltest)
  * Phylogeny (ML or Bayesian)
  * ML tree with bootstrap values 
  * Bayesian tree has posterior values

* **Bootstrapping:** Bootstrapping is a resampling analysis that involves taking columns of characters out of your analysis, rebuilding the tree, and testing if the same nodes are recovered. This is done through many (100 or 1000, quite often) iterations. If, for example, you recover the same node through 95 of 100 iterations of taking out one character and resampling your tree, then you have a good idea that the node is well supported (your bootstrap value in that case would be 0.95 or 95%). If you get low support, that suggests that only a few characters support that node, as removing characters at random from your matrix leads to a different reconstruction of that node.

* **Bayesian posterior probabilities:** A bayesian tree has posterior values but if we can also add bootstrap values for each branch along with posterior values from the newick tree.

