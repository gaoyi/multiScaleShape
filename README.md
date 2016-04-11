# README #

### Summary ###
Following the paper:

Gao, Yi, Benjamin Corn, Dan Schifter, and Allen
Tannenbaum. "Multiscale 3D shape representation and segmentation with
applications to hippocampal/caudate extraction from brain MRI."
Medical image analysis 16, no. 2 (2012): 374-385.

this repository contains code to perform multi-scale shape
representation and decomposition using 3D wavelet transformation.

### What is this repository for? ###

* Forward and inverse Wavelet transformation of 2D and 3D ITK images.
* Shape decomposition, reconstruction, and filtering. This includes hard and soft thresholding

### How do I get set up? ###

* CMake (cmake.org)
* ITK (itk.org, Compilation needed)
* GNU gsl library (http://www.gnu.org/software/gsl, Compilation needed)
* boost library (download from www.boost.org, un-zip, NO compilation needed)

### Who do I talk to? ###

* Yi Gao (gaoyi@gatech.edu)

### Usage agreement ###

In scientific publications (journal publications, conference papers,
technical reports, presentations at conferences and meetings) that use
this code, you should cite the following paper:

Gao, Yi, Benjamin Corn, Dan Schifter, and Allen
Tannenbaum. "Multiscale 3D shape representation and segmentation with
applications to hippocampal/caudate extraction from brain MRI."
Medical image analysis 16, no. 2 (2012): 374-385.


### To do ###

* Due to historical reason, there is some un-necessary data conversion through cArray?D. This should be removed.
* The dependence on the boost library is only the shared pointer. It can also be removed.
