AdelaideRMF: Robust Model Fitting Data Set

AdelaideRMF is a data set for robust geometric model fitting (homography estimation and fundamental matrix estimation). We collected a set of image pairs and manually labelled the keypoint correspondences which were obtained by SIFT matching.

1. adelaidermf folder contains Matlab .mat files.
Each mat file contains the following:	
(a) img1 - left image
(b) img2 - right image
(c) data - keypoint correspondences (x1,y1,1,x2,y2,1) where (x1,y1) in img1 and (x2,y2) in img2
(d) score - SIFT correspondences matching score
(e) label - 0 indicates gross outliers and others indicate the structure membership

2. adelaidermf_images folder contains original camera images.

[Homography Data]
barrsmith, johnsona, oldclassicswing, johnsonb, physics, ladysymon, sene, elderhalla, library, elderhallb, napiera, unihouse, bonhall, napierb, unionhouse, bonython, neem, hartley, nese

[Fundamental Matrix Data]
breadcartoychips, cubechips, biscuit, breadcube, cubetoy, biscuitbook, breadcubechips, dinobooks, biscuitbookbox, breadtoy, toycubecar, boardgame, breadtoycar, carchipscube, game, cube, gamebiscuit, book, cubebreadtoychips


To use this dataset, please cite this paper:  
"Hoi Sim Wong, Tat-Jun Chin, Jin Yu, and David Suter, Dynamic and Hierarchical Multi-Structure Geometric Model Fitting, International Conference on Computer Vision (ICCV), Barcelona, Spain, 2011."

For enquires regarding this data set, please email hoisimwong@gmail.com

