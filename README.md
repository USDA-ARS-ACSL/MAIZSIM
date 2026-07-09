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

Updates since 2025
Date – Comment
7/7/2026 12:17:33 PM – increased max fert times to 175
3/23/2026 3:11:46 PM – updated version number in splash screens
3/20/2026 9:52:28 AM – redoing NH4 related changes - used the wrong repo for source previously
3/20/2026 9:23:31 AM – fixed issue where NH4 was set to 0 after input because BD had not been read in yet
9/24/2025 11:59:52 AM – added a comment
9/24/2025 11:59:14 AM – harmonized the single/double precision calculations. was getting NaN when dt was declared a double precision. Declared more variables and removed some that were no longer needed
8/12/2025 5:09:25 PM – added code to catch case where there was too much water to infiltrate in one day
7/23/2025 4:10:22 PM – removed or modified code to compile in linux and removed outdated fortran expressions
5/13/2025 9:53:02 AM – added date to header
5/13/2025 9:52:26 AM – added double descriptor to call for abio to allow it to work in linux
5/13/2025 9:51:05 AM – renamed weather structure to make it compatible with linux (hourly weather had this same named block)
5/13/2025 9:33:46 AM – changed weather common block name to be compatible with linux
3/7/2025 9:13:15 AM – added more comments
1/15/2025 1:18:14 PM – when temps went <0 only minimum temperature was constrained. Added constraint to 2C for max temp also. Also constrained TminT
