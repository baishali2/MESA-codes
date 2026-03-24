#MESA-codes:
This is my MSC project repository, which has 3 directories added: the original work directory of MESA, and 2 other directories which have been used to run CE evolutions in MESA.

mesa-24.08.1/star/work : the original work directory of MESA version-24.08.1, which has been used in this project. Can be used for evolution of a star, masses, compositions and various other properties can be varied.

work_dir_abha : Specific CE evolution code developed by my senior Abha. Uses a pre-built model of a 12 Msun RGB star in helium burning phase. To use this code, the user has to copy the run_star_extras.f90, and inlist_project file, replace the default files from the work directory, and run MESA.

work_dir_baishali : Contains my contribution to the CE evolution code. The src folder contains 2 run_star_extras, one for my RGB run and one for my AGB current run, the AGB version contains code for inputting drag force data from 3D simulation. The RGB run uses a 2 Msun star evolved upto RGB phase, and the AGB run uses the same 2 Msun star evolved upto AGB phase. The mass of the companion can be varied in the run_star_extras.f90.
