#MESA-codes:
This is my project repository, which has 3 directories added: the original work directory of MESA, and 2 other directories which have been used to run CE evolutions in MESA.

mesa-24.08.1/star/work : the original work directory of the version of MESA used in this project. Can be used to simulate a star of various masses, compositions and other desired settings.

work_dir_abha : Specific CE evolution code developed by my senior Abha. Uses a pre-built model of a 12 Msun RGB star in helium burning phase. To use this code, the user has to copy the run_star_extras.f90, and inlist_project file, and run MESA with the default files replaced with those.

work_dir_baishali : Contains my contribution to the CE evolution code. The src folder contains 2 run_star_extras, one for my RGB run and one for my AGB current run, the AGB version contains code for inputting drag force data from 3D simulation.
