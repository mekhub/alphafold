# AlphaFold
(C) R. Das, Stanford University 2018

## What this is
A Python package for modeling the statistical physics of RNA folding at the secondary structure level. 

Goals:
 * Code that is easy to read so that humans can easily _extend_ it to model new RNA physics 
 * Code with numerous tests built in so that extensions are _correct_
 * A package that can learn from the huge data sets our lab is collecting.

A separate C++ package `alphafoldplus` with the same functionality and matching python bindings and likely up to 100x the speed is being developed separately in a [private repository](https://github.com/rhiju/alphafoldplus).

## Features 
This code brings together features pioneered in (but scattered across) prior packages:
 * Multi-strand calculations
 * Circular RNAs
 * Co-axial stacking
 * True partition function calculations in `N^3` time
 * Base pair probability estimates
 * Gradients of predicted observables with respect to energy model parameters, to enable learning from data
 * Enumerative backtracking to get all structures and their Boltzmann weights
 * Stochastic backtracking to get Boltzmann-sampled structures
 * Minimum free energy structures
 * Generalized base pairs (e.g., both Watson-Crick and Sugar/Hoogsteen G-A pairs) (_coming soon_)
 * 'Classic' Turner2004 & ContraFold parameters (_coming soon_)
 * Modeling of ligand/protein binding to RNA hairpins and internal loops (_coming soon_)
 * Modeling of protein binding to RNA single-stranded segments (_coming soon_)
 
This code also presents entirely new features, based on recent theoretical insights from R. Das & students:
 * Cross-checks based on computation of the partition function `N` different ways for each RNA.
 * Loop penalties that rise like the logarithm of the number of loop nucleotides, still in `N^3` time (_coming soon_)
 * Parameters for chemically modified bases, and some modified backbones, based on Rosetta calculations (_coming soon_)
 * Linear motifs identified by Rosetta or by crystallography as having favorable energy bonuses (_coming soon_)
 * Modeling of protein binding to RNA, including proper steric exclusion effects. (_coming soon_)
 * Modeling of RNA tertiary contacts, through a novel iterative sampling method, Rosetta-calculated properties of the contacts, and efficient C_eff calculations. (_coming soon_)
 * Tracking and propagation of estimated model uncertainties. (_coming soon_)
 * Efficient learning from large data sets through stochastic gradient-based calculations (_coming soon_)
 * Easy install through `sudo pip` (_coming soon_)
 
## License
This code is being released with the MIT license. So you can distribute it with your code. 

## Getting started
Clone this repository, and just type:
```
./alphafold.py
```
to run tests on a bunch of example sequences.

To run on tRNA(phe) from yeast:
```
./alphafold.py -s GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA
```

To circularize:

``` 
./alphafold.py -s GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA --circle
```

To run on a multi-strand system, type:
```
./alphafold.py -s GCAACG CGAAGC
```

To re-run tRNA in a totally weird way:
```
./alphafold.py -s UGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGAC --circle
```
Should get the same answer as above linear case!

## Contributing
More information on making contributions coming soon.
