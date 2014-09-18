MAIZSIM is a corn simulation model developed by the USDA-ARS Crop Systems and Global Change Laboratory and Univ. of Washington School of Environmental and Forest Sciences. 

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

More documention is being prepared. 

Input files for one test case are in the folder Test_Ver9.1


