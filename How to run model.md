Maizsim has to be run from the command line. A "run" file will contain paths and filenames of the input files:
D:\MAIZSIM07\MDEasternShore\Del06\DEL06.wea
D:\MAIZSIM07\MDEasternShore\Del06\DEL06.tim
D:\MAIZSIM07\MDEasternShore\Del06\BiologyDefault.bio
D:\MAIZSIM07\MDEasternShore\Del06\WyeClimate.dat
D:\MAIZSIM07\MDEasternShore\Del06\Del.nit
D:\MAIZSIM07\MDEasternShore\Del06\NitrogenDefault.sol
D:\MAIZSIM07\MDEasternShore\Del06\WyeSoil.soi
.
.
.
In the examples, the run file is always prefixed with 'run' and has dat  as an extension, i.e., "RunDel06.dat". The model is run on the command from a DOS prompt. Assuming your executable is in the folder d:\maizsim07\MDEasternShore and the input data are in subfolders - d:\maizsim07\MDEasternShore\DEL06 and d:\maizsim07\MDEasternShore\DEL07 for example, you would run the model as:

D:\Maizsim07\MDEasternShore>2dsoil .\DEL06\runDEL06.dat

the .\ is a relative path address and will tell the operating system to look for the run file in a subdirectory called wye06. The model executable (exe file) is still called 2dsoil because of how the compiler was originally set up. We plan to change it in the future. The folder with the 2dsoil.exe file should also have the files crop.dll and lightenv.dll. There are 5 output files:

DEL06.g01  -- plant output
DEL06.g02  --detaile leaf output
DEL06.G03  --water, temperature, concentration values at the nodes of the soil grid
DEL06.G04  --root information for the nodes or elements
DEL06.G05  -- surface fluxes, et, rain, transpiration
DEL06.G06  bottom and top boundary fluxes

I usually make a folder for each simulation and keep all the related files in that folder. One can put some files in other folders to reduce duplication, for example put all the grid files into a grid folder. Initially, I used to put all the related data in the same folder, for example the weather files in a weather folder or the variety files in a variety folder. Since the full file paths are specified in the run.dat file, the files can be located anywhere. 
