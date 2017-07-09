
# The Random Cluster Model for Robust Geometric Fitting

This package contains the source code which implements robust geometric model fitting proposed in:

T.T. Pham, T.-J. Chin, J. Yu and D. Suter
The Random Cluster Model for Geometric Model Fitting  
In Proc. IEEE Conf. on Computer Vision and Pattern Recognition (CVPR),  Providence, Rhode Island, USA, 2012.

T. T. Pham, T.-J. Chin, J. Yu and D. Suter
The Random Cluster Model for Robust Geometric Fitting
IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI).

Copyright (c) 2012 Trung T. Pham and Tat-Jun Chin
School of Computer Science, The University of Adelaide, South Australia
The Australian Center for Visual Technologies
http://www.cs.adelaide.edu.au/~{trung,tjchin}

If you encounter any issues with the code, please feel free to contact me at:
trung.pham@adelaide.edu.au

Last updated: 09 July 2017.

————————-
0. Libraries
—————————
This software uses the Multi-label optimization toolbox developed by Olga Veksler and Andrew Delong, which can be downloaded from http://vision.csd.uwo.ca/code/gco-v3.0.zip. We include this toolbox to our package.

This program also makes use of Peter Kovesi and Andrew Zisserman's MATLAB functions for multi-view geometry
(http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/ http://www.robots.ox.ac.uk/~vgg/hzbook/code/).
 
----------------------------
1. Installation Instructions
----------------------------
* Uncompress the package.
* Install GCO library. 
..* Go to gco-v3.0/matlab directory. 
..* Run GCO_UnitTest.m. The mex file should be compiled automatically. For more information, please see readme.txt file under gco-v3.0/matlab directory. 
* Run 'make.m' file.

---------------
2. Run examples
---------------
* Run homo_est.m for multiple planar homography detection.
* Run funda_est.m for multiple 2-view motion segmentation.

Users are able to test other datasets by changing the data name manually in the funda_test.m or homo_test.m files. 

Note: We have tested the code under Linux (Ubuntu) and Mac.  


-----------
3. License
-----------
The program is free for non-commercial academic use. Any commercial use is strictly prohibited without the authors' consent. Please acknowledge the authors by citing the above paper in any academic publications that have made use of this package or part of it.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


