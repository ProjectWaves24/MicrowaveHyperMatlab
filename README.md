# MicrowaveHyperMatlab

 Matlab code (located in  zip-files) is used for real-time temperature monitoring in the process
 of microwave thermometry.
 
1) The code in zip file MatlabCodeAllTimes.zip is
used for generation of data for C++/PETSc code located in  private directory

https://github.com/ProjectWaves24/MESH

This  code was used for computations of the paper

L. Beilina, M. G. Aram, E. M. Karchevskii,
An adaptive finite element method for solving 3D electromagnetic volume integral equation with applications in microwave thermometry,
Journal of Computational Physics, V. 459, 2022, 111122, ISSN 0021-9991,
https://doi.org/10.1016/j.jcp.2022.111122.


as well as to obtain least squares reconstructions reported in the paper

 M. G. Aram, L. Beilina, H. Dobsicek Trefna, Microwave Thermometry with Potential Application in Non-invasive Monitoring of Hyperthermia, Journal of Inverse and Ill-posed problems, 2020
DOI: https://doi.org/10.1515/jiip-2020-0102

Read file README_MATLAB for further details about the code.  Note that C++/PETSc code for adaptive reconstruction is located in private directories: to be able get access write enquiry to larisa@chalmers.se

2) The code located in MATLAB_TSVD.zip or in the directory MATLAB_TSVD/
uses truncated SVD in order to obtain reconstructions. This code is used for all computations of the paper

L. Beilina, E.Lindström, L. Frischauf, Daniel Mc Kelvey, Truncated SVD for applications in microwave thermometry, to appear in Springer Proceeding in Mathematics and Statistics.
