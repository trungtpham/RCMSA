
# The Random Cluster Model for Robust Geometric Fitting

This package contains the source code which implements robust geometric model fitting proposed in:

T.T. Pham, T.-J. Chin, J. Yu and D. Suter
The Random Cluster Model for Geometric Model Fitting  
In Proc. IEEE Conf. on Computer Vision and Pattern Recognition (CVPR),  Providence, Rhode Island, USA, 2012.

T. T. Pham, T.-J. Chin, J. Yu and D. Suter
The Random Cluster Model for Robust Geometric Fitting
IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI).

Related papers:

T. T. Pham, T.-J. Chin, K. Schindler and D. Suter,
Interacting Geometric Priors for Robust Multi-Model Fitting
IEEE Transactions on Image Processing

T. T. Pham, T.-J. Chin, J. Yu and D. Suter
Simultaneous Sampling and Multi-Structure Fitting with Adaptive Reversible Jump MCMC
In NIPS 2011, Granada, Spain.

Copyright (c) 2012 Trung T. Pham and Tat-Jun Chin
School of Computer Science, The University of Adelaide, South Australia
The Australian Center for Visual Technologies
http://www.cs.adelaide.edu.au/~{trung,tjchin}

If you encounter any issues with the code, please feel free to contact me at:
trung.pham@adelaide.edu.au

Last updated: 15 Dec 2017.

----------------------------
0. Libraries
----------------------------
This software uses the Multi-label optimization toolbox developed by Olga Veksler and Andrew Delong, which can be downloaded from http://vision.csd.uwo.ca/code/gco-v3.0.zip. We include this toolbox to our package.

This program also makes use of Peter Kovesi and Andrew Zisserman's MATLAB functions for multi-view geometry
(http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/ http://www.robots.ox.ac.uk/~vgg/hzbook/code/).
 
----------------------------
1. Installation Instructions
----------------------------
* Uncompress the package.
* Install GCO library. 
	- Go to gco-v3.0/matlab directory. 
	- Run GCO_UnitTest.m. The mex file should be compiled automatically. For more information, please see readme.txt file under gco-v3.0/matlab directory. 
* Run make.m file.

---------------
2. Run examples
---------------
* Run homo_eval.m to test the method using AdelaideRMF dataset.
* Run funda_eval.m to test the method using AdelaideRMF dataset.

Note: We have tested the code under Ubuntu 16.04 and Matlab R2017a.


