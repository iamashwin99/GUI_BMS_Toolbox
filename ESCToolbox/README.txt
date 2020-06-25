This archive includes files necessary to create and work with enhanced-self-
correcting cell models in MATLAB. The archive is organized into two principal
folders: OCV_FILES and DYN_FILES.

Inside OCV_FILES, you will find data collected from six different cells that
can be used to create an open-circuit-voltage relationship for these cells.
The data are contained in the subfolders A123_OCV, ATL_OCV, and so forth.
The "raw" data from the cell tester is stored in an Excel file format; it is
also converted into a more convenient ".mat" MATLAB file format in each of the
subfolders. In the main folder, you will find the following ".m" files: 
- makeMATfiles.m:  Converts the Excel files into .mat files
- plotMATfiles.m:  Plots the different OCV script voltages for the .mat files
- processOCV.m:    Creates an OCV relationship for one cell
- runProcessOCV.m: Creates OCV relationships for all cells.
You should be able to execute runProcessOCV.m directly to create the OCV 
relationships, which will be stored as .mat files (A123model-ocv.mat, etc.)

Inside DYN_FILES, you will find data collected from the same six cells that 
can be used to create a dynamic relationship for these cells.
The data are contained in the subfolders A123_DYN, ATL_DYN, and so forth.
The "raw" data from the cell tester is again stored in Excel file format; it 
is also converted into a more convenient ".mat" MATLAB file format in each of
the subfolders. In the main folder, you will find the following ".m" files:
- OCVfromSOCtemp.m: A utility function to compute OCV from SOC for final model
- SOCfromOCVtemp.m: A utility function to compute SOC from OCV for final model
- getParamESC.m:    A utility function to get any desired model parameter value
- makeMATfiles.m:   Converts the Excel files into .mat files
- plotMATfiles.m:   Plots the different dynamic script voltages from .mat files
- processDynamic.m: Creates the dynamic relationship for one cell
- runProcessDynamic.m: Creates the dynamic relationship for all cells
- setupDynData.m:   A utility script to define common variables
- simCell.m:        A utility function to simulate a cell's response to stimulus
You should be able to execute runProcessDynamic.m directly to create the 
dynamic relationships (after having run the runProcessOCV.m script first).
These will be stored as .mat files (A123model.mat, etc.)

---

All files in this archive are copyright (c) 2015 by Gregory L. Plett of the 
University of Colorado Colorado Springs (UCCS). This work is licensed under 
a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Intl. License, 
v. 1.0. It is provided "as is", without express or implied warranty, for 
educational and informational purposes only.

These files are provided as a supplement to: Plett, Gregory L., "Battery
Management Systems, Volume I, Battery Modeling," Artech House, 2015.
