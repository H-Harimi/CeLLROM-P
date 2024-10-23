This archive includes files necessary to create and work with reduced-order 
models of lithium-ion battery cells. The archive is organized into three principal
categories: 1) creating ROMs, 2) simulating ROMs, and 3) simulating FOMs.

1) The files associated with creating ROMs are described here. 

The parameter values describing the cell to be modeled with a ROM are stored in 
an Excel spreadsheet. “Doyle_parameter_list.xlsx” is an example of the correct 
format for  the “Doyle” cell described in Doyle, Fuller, and Newman’s 1996 paper. 
The DRA is used to create a ROM, and the DRA tuning parameters are also stored
in an Excel spreadsheet.  “Doyle_DRA_list.xlsx” is an example of the correct
format that can be used with the Doyle cell.

To create a ROM, use the runGenROM.m function. It shows an example of how to
call the “genROM.m” code, which does the work.  This, in turn, calls the “dra.m”
function, which is the core DRA procedure.

The ROM can contain as many or as few transfer functions as you like, evaluated
at whatever points in the cell as you like. (To be able to compute cell voltage,
there are some required transfer functions, but genROM does not check for this, 
so you can make quite general ROMs.) You can implement your own transfer functions
using the provided ones as a template. The key requirements are that the transfer-
function “.m” files have the same output structure — the input structure is more 
flexible as it is specified by the user in the DRA tuning-parameter Excel file.  
The example transfer functions provided in this archive are:
- tf_ce.m:        Electrolyte concentration (less ce0)
- tf_cse.m:       Solid surface concentration (less cs0)
- tf_gradce.m:    Gradient of electrolyte concentration
- tf_gradphie1.m: Gradient of linear part of electrolyte potential
- tf_gradphis.m:  Gradient of solid potential
- tf_j.m:         Reaction flux
- tf_phie1.m:     Linear part of electrolyte potential
- tf_phis.m:      Solid potential (less v(t) in positive electrode)
- tf_phise.m      Interfacial potential difference

The other related functions are:
- dra.m:            The core DRA code itself
- genROM.m:         The function that creates the ROMs
- readDRATable.m:   A utility procedure for reading the Excel DRA tuning params
- readParamTable.m: A utility procedure for reading the Excel parameter table
- runGenROM.m:      Example of how to use “genROM.m”

Note that the DRA can be quite memory and computationally intensive. Some pre-
computed ROMs are provided in the folders “Doyle_ROM_Files_4” for an example 
ROM having 4 general-purpose states plus one integrator state, and in
“Doyle_ROM_Files_9” for an example ROM having 9 general-purpose states plus one
integrator state. The corresponding Excel parameter and DRA files are included
in these folders so that you can see how they were created (however, they require
about 64GB of physical RAM plus a lot of patience).  The example in the parameter
and DRA Excel spreadsheets in the main level folder of the archive can be run
with much more modest computational resources, although the results will not be
quite as good.

2) The files associated with simulating ROMs are
- runSimROM.m: This is an example of how to use “simROM.m” for a number of 
               example cases
- simROM.m:    This simulates a reduced-order model for some input stimulus

3) The files associated with simulating FOMs are
- runSimFOM.m: This is an example of how to use “simFOM.m” for the same example
               cases as for “runSimROM.m”
- simFOM.m:    This simulates the full-order model of a cell, using cell params
               from an Excel spreadsheet (same as for genROM.m). Note that this
               code requires COMSOL software with the MATLAB LiveLink toolbox
               in order to run.

Since not everyone has the COMSOL software, example outputs from simFOM are 
stored in the folder “FOM_Output_Files.”

Finally, a simple graphical user interface, compareROMvsFOM.m (and .fig), is 
provided to be able to visualize results from simROM and simFOM.

---

All files in this archive are copyright (c) 2015 by Gregory L. Plett of the 
University of Colorado Colorado Springs (UCCS). This work is licensed under 
a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Intl. License, 
v. 1.0. It is provided "as is", without express or implied warranty, for 
educational and informational purposes only.

These files are provided as a supplement to: Plett, Gregory L., "Battery
Management Systems, Volume I, Battery Modeling," Artech House, 2015.
