![alt tag](https://travis-ci.org/daviddesancho/RepeatDesigner.svg?branch=master)
![alt tag](https://landscape.io/github/daviddesancho/RepeatDesigner/master/landscape.svg?style=flat)
[![Coverage Status](https://coveralls.io/repos/github/daviddesancho/RepeatDesigner/badge.svg?branch=master)](https://coveralls.io/github/daviddesancho/RepeatDesigner?branch=master)


# RepeatDesigner
A library for re-designing repeat proteins

Functionality
-------------
The package uses homology modelling tools for re-designing repeat proteins.
Two different scenarios are envisioned.

* Concerted engineering of repeat residues, so that mutations are introduced
simultaneously in all repeats.
* Decoupled engineering of repeats, in order to produce different specificities.


Dependencies
------------

* BioPython - http://biopython.org/
* Modeller - http://salilab.org/modeller/
* Seaborn - https://github.com/mwaskom/seaborn

These dependencies will be checked by the installation process with the
 exception of Modeller, which requires prior installation by the end-user.
