#MESA-codes:
This is my MSC project repository, which has 3 branches: the *main* branch, the *work_dir_abha* branch, the *work_dir_baishali* branch, the last two contain directories of the same name.

*main* branch : contains mesa-24.08.1/star/work : the original work directory of MESA version-24.08.1, which has been used in this project. Can be used for evolution of a star, masses, compositions and various other properties can be varied.

*work_dir_abha* branch : Contains directory of the same name. Specific CE evolution code developed by my senior Abha. Uses a pre-built model of a 12 Msun RGB star in helium burning phase. To use this code, the user has to copy the run_star_extras.f90, and inlist_project file, replace the default files from the work directory, and run MESA.

*work_dir_baishali* branch : Contains my contribution to the CE evolution code. 2 versions of inlist_project for RGB and AGB simulation are added. src folder contains all models of the AGB simulation code, the details of these models is added into the description of the directory.
