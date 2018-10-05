This is the read me file for the Wrapper (/home/user1/Wrapper) this directory contains the list of files and folders that will that are conatained within this directory and their function.

IMPORTANT NOTE: when entering the pmt number only input the numbers that are different e.g. NB0081 should be inputted as 81 and NB0166 should be inputted as 166

IMPORTANT NOTE: To run Geco enter sudo su into the terminal, then type /usr/local/bin/CAENGECO2020 (su /usr/local/bin/CAENGECO2020 should also work) 

in Wrapper: 
HVScan.txt -> A text file that contain the ids and voltages for each pmtm the first column is the id. 5 voltages are listed in 100 volt bins, the 5 voltages are a spread around the the voltage in the final column. The final entry is value hamamatsu determined for each PMT. 

Setup.sh -> setup for root (unsure why this is here)

wave_fitter .cpp .o and .d files -> Main, object and dictory files for the wave_fitter program, the main is the only one that is edited, the other two are for the computer to build the program sucessfully. 

wave_fitter -> Unknown complied program, makes histograms continiously 

wave_process -> Unknown complied program, makes 3 empty histograms appear

spe2.root -> Uknown file, empty, created from wave_process

And the folders for SPE_Analysis, Dark_rate, After_Pulsing, After_Pulsing_PICO and Gain_Test

######### new directory #########

in Wrapper/SPE_Analysis:

SPE_Gen and SPE_Analysis .cpp .o and .d files -> Main, object and dictory files for the SPE_Gen and SPE_Analysis programs, the main is the only one that is edited, the other two are for the computer to build the program sucessfully. 

SPE_Gen -> Complied program for simple singe photo peak histogram to be created, enter in pmt number and voltage when prompted

SPE_Analysis -> Unknown complied program, gives peak to valley ratio and area and centroid.  Enter in pmt number and voltage when prompted

SPE_root -> results folder for SPE_Gen and SPE_Analysis

######### new directory #########

in Wrapper/Dark_Rate:

SPE_Gen and SPE_Fit and wave_fitter .cpp .o and .d files -> Main, object and dictory files for the SPE_Gen and SPE_Fit and wave_fitter programs, the main is the only one that is edited, the other two are for the computer to build the program sucessfully.

SPE_Gen ->  Complied program for simple singe photo peak histogram to be created, but takes much longer than the one in Wrapper/SPE_Analysis. Enter in pmt number and voltage when prompted

SPE_Fit -> Unknown complied program, fits the dark noise taken by the SPE_Gen program, however cannot seem to fit the data, outputs the Dark Counts.  Enter in pmt number and voltage when prompted.

wave_fitter -> Unknown complied program, produces histograms of waves constantly 

SPE -> results folder for the SPE_Gen and SPE_Fit

######### new directory #########

in Wrapper/Gain_Test:

SPE_Gen and SPE_Fit and SPE_Roll .cpp  files ->  Main files for the SPE_Gen and SPE_Fit and SPE_Roll, object and dictionary files not present for some reason. 

SPE_Gen -> Complied program for simple singe photo peak histogram to be created, enter in the PMT number and then the relevent coloumn from HVScan.txt to enter in the voltages 

SPE_Fit -> Unknown complied program, crashes with segmentation violation when attempted to be used. Enter in the PMT number and then the relevent coloumn from HVScan.txt to enter in the voltages 

SPE_ROLL -> Uknown complied program

HV_SPE -> Results Folder, defiantly for the SPE_Gen, probably for the SPE_Fit, SPE_ROLL too once they work 

######### new directory #########

in Wrapper/After_pulsing:

Unknown directory, assumed to be unfinished at time of writing 18/09/2018

######### new directory #########

in Wrapper/After_pulsing_PICO:

wave_fitter.cpp  wave_process.cpp ->  Main files for the wave_fitter  wave_process programs object and dictionary files not present for some reason.

wave_fitter-> complied program, produces continious spectra, presumably from the PICO

wave_process-> complied program, produces a continious stream of enteries, presumably fomr the PICO 

setup.sh -> setup for root 

SPE -> unknown folder, for results?
