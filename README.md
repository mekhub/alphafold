# AlphaFold
(C) R. Das, Stanford University 2018

## What this is
A Python package for modeling the statistical physics of RNA folding at the secondary structure level. 

Goals:
 * Code that is easy to read so that humans can easily _extend_ it to model new RNA physics 
 * Code with numerous tests built in so that extensions are _correct_
 * A package that can learn from the huge data sets our lab is collecting.

A separate C++ package `alphafoldplus` with the same functionality and matching python bindings and likely up to 100x the speed is being developed separately in a private repository.

## Features 
This code brings together features pioneered in (but scattered across) prior packages:
 * Multi-strand calculations
 * Circular RNAs
 * Co-axial stacking
 * True partition function calculations in `N^3` time
 * Base pair probability estimates
 * Gradients of predicted observables with respect to energy model parameters, to enable learning from data
 * Enumerative backtracking to get all structures and their Boltzmann weights (_coming soon_)
 * Stochastic backtracking to get Boltzmann-sampled structures (_coming soon_)
 * Minimum free energy structures (_coming soon_)
 * 'Classic' Turner2004 & ContraFold parameters (_coming soon_)
 * Generalized base pairs (e.g., both Watson-Crick and Sugar/Hoogsteen G-A pairs) (_coming soon_)
 
This code will also present entirely new features, based on recent theoretical insights from R. Das & students:
 * Cross-checks based on computation of the partition function `N` different ways for each RNA.
 * Loop penalties that rise like the logarithm of the number of loop nucleotides, still in `N^3` time
 * Parameters for chemically modified bases, and some modified backbones, based on Rosetta calculations
 * Linear motifs identified by Rosetta or by crystallography as having favorable energy bonuses
 * Modeling of RNA tertiary contacts, through a novel stochastic sampling method and Rosetta-calculated properties of the contacts.
 * Tracking and propagation of estimated model uncertainties.
 * Efficient learning from large data sets through stochastic gradient-based calculations
 * Easy install through `sudo pip`
 
## License
This code is being released with the MIT license. So you can distribute it with your code. 
