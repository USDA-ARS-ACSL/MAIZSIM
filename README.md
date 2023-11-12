MAIZSIM is a corn simulation model developed by the USDA-ARS Adaptive Cropping Systems Laboratory and Univ. of Washington School of Environmental and Forest Sciences. Note that the dev branch is used to hold recent code. After testing, it will be merged with the master branch.  If you want to collaborate, please fork the code, make your changes and  then create a pull request.

The developers are 

Soo-Hyung Kim of the Univ of Washington
Dennis Timlin, David Fleisher, and V.R. Reddy of the USDA-ARS

others who have collaborated include:
Yang Yang, now at Dow Agrosciences

Annette Dathe Cornell University, Ithaca, NY
 diffusive root model
 
Jong-Ahn Chun APEC Climate Center, Korea
 CO2 and water
 
Sahila Beegum, Univ of Nebraska
 gas transport and respiration
 
Wenguang Sun, Colorado State Univ.
 gas transport and respiration
 

MAIZSIM is a mechanistic model of maize growth, development and yield. It is written in C++ (crop) and FORTRAN (soil). 

The model is interfaced with 2DSOIL, a two dimensional simulator of soil water and heat movement, and solute transport. This model is written in FORTRAN and is the main model. 2DSOIL calls the crop model as a subroutine. There are two subprojects, Crop Source and Soil Source.

The code compiles in visual studio.net. We used Intel Fortran (OneAPI-2023.2) and Visual Studio Professional 2022. Macros will copy the compiled libraries dll's for the crop mode to the folder with the soil's code if you make sure to keep your soil and crop source code in the original folders ('crop source' and 'soil source') which are the two subprojects.

It has recently been revised to compile under Linux by Kyungdahm Yun of the Univ. of Washington. See below for a link to a docker image with a makefile and the compilers. There is a makefile and instructions contained in a docker image. The link is at this repository. The link contains a docker image so you can test it in a Windows environment using WSL2.

https://github.com/precision-sustainable-ag/BuildMaizsim

More documention is being prepared. See the "how to run model" file for information on how to set up the input files and run the executable from the command line. 


See the ExcelInterface repository for an excel based interface and example input files.
https://github.com/USDA-ARS-ACSL/ExcelInterface

The most recent updates for the excel interface work with the most recent version of maizsim
test for readme. Use the same tag number



