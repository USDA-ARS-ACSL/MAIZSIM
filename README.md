MAIZSIM is a corn simulation model developed by the USDA-ARS Adaptive Cropping Systems Laboratory and Univ. of Washington School of Environmental and Forest Sciences. 

The developers are 

Soo-Hyung Kim of the Univ of Washington
Dennis Timlin, David Fleisher, and V.R. Reddy of the USDA-ARS

others who have collaborated include:
Yang Yang, now at Dow Agrosciences
Annette Dathe Norwegian University of Life Sciences
Jong-Ahn Chun APEC Climate Center, Korea

MAIZSIM is a mechanistic model of maize growth, development and yield. It is written in C++ (crop) and FORTRAN (soil). 

The model is interfaced with 2DSOIL, a two dimensional simulator of soil water and heat movement, and solute transport. This model is written in FORTRAN and is the main model. 2DSOIL calls the crop model as a subroutine. 

The code compiles in visual studio.net. It has recently been revised to compile under Linux by Kyungdahm Yun of the Univ. of Washington 

There are two subprojects, Crop Source and Soil Source. 

More documention is being prepared. See the "how to run model" file for information on how to set up the input files and run the executable from the command line. 

Input files for test cases related to the data from the Maryland Eastern Shore and Delaware used in the first paper are in the folder Test_Ver9.1 and MDEasternShore. The MDEasternShore data work with the newest version (1.42) of the crop model as of Oct 25, 2016. These files have not been updated for the new version (yet).

The SolarCorr.Zip file contains intput and output corresponding to the 9/6/2018 update. I still have to fix the other input files so they work with this version.


