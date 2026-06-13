To run any of the simulations: There are two ways to run the simlations. First, you need to have the original mesa-24.08.1 or any other latest version of MESA on your system/server. The only modifications over original MESA work directory are contained in these 2 files: inlist_project and run_star_extras. 

Then, either you can download this entire directory and remove the versions of inlist_project and run_star_extras that you would not need. For example, to run the RGB simulation, you will need the starting model rgb_model_558_withbc.mod, inlist_project_rgb, and run_star_extras_rgb.f90 from the src folder. You can remove the inlist_project_agb, and all the versions of run_star_extras_agb. **Don't forget to change the names of these files to simply inlist_project and run_star_extras.f90, you don't need to change the starting model's name**. Then you can run it following MESA's standard procedure. 

Or, following MESA's standard procedure, you can copy the mesa-24.08.1/star/work directory somewhere outside the mesa directory.Then, you can download the versions of inlist_project, run_star_extras and the starting model that you will need, and replace the original MESA versions of inlist_project and run_star_extras with these versions. Then again you can run following MESA's standard procedure. For example, you can just download inlist_project_rgb, run_star_extras_rgb.f90, and rgb_model_558_withbc.mod, and delete the original versions of inlist_project and run_star_extras.f90. **In this case too, don't forget to change the names of these files to simply inlist_project and run_star_extras.f90, you don't need to change the starting model's name.**

Added mention: you cannot have more than 1 run_star_extras in the src folder of the work directory. So when you need 1 run_star_extras, remove all the others.

Details of the different files here:


inlist_project_rgb : inlist project for the RGB simulation.

inlist_project_agb : inlist project for the AGB simulation.

rgb_model_558_withbc.mod : 2 Msun star evolved up to RGB phase, used in the RGB simulation. This model was made with the boundary conditions on - as instructed by Abha.

agb_without_bc_9765.mod : 2 Msun star evolved upto AGB phase, used in the AGB simulation. This model was made without the boundary conditions, it was seen that no boundary conditions numerically stabilises the simulation, otheerwise the simulation crashes very early.

src: This folder contains all the run_star_extras for different models and different simulations. The run_star_extras file contains all of the physics related to the CE phase, which have been input as user-defined routines.

The src folder contains:

1. run_star_extras_rgb.f90 : The main code for the RGB final simulation.

2. run_star_extras_agb_model1.f90 : Code for AGB simulation Model 1. I calculated the drag force exerted by the envelope gas on the companion in the frame of the core particle. After that, I
calculated the ϕ component of this force(Fϕ). This Fϕ
is input to the simulation as drag force. The force data is read as a function
of time and interpolated to find the drag force in each MESA timestep.

3. run_star_extras_agb_model2.f90 : Code for AGB simulation Model 2. This model has the same input drag force as that of model 1, but with a time shift of 20.3 days. Since I placed the companion inside the envelope at 122 R⊙ and
not at 124 R⊙ as in the 3D simulation, so reading the drag force data from the
beginning introduces a time mismatch, relative to the position of the companion.
To account for this, I determined the time in the 3D simulation when the orbital
separation reached 122 R⊙, which was 20.3 days. I then started reading the drag
force data from 20.3 days onward, keeping the start time of my simulation at
zero.

4. run_star_extras_agb_model3_final.f90 : Code for AGB simulation Model 3. In the first part of the simulation, up to 97 days, I neglect the motion of the
core with respect to the envelope. Consequently, by ignoring the motion of the core, I only calculate the ϕ compo-
nent of the drag force on the companion(F2,ϕ). From 97 days onwards, I again take into account the motion of the core.
Accordingly, I take the drag force to be Fϕ, calculated with the same method
as in model 1. To solve the time mismatch issue, instead of reading the drag force data from
a later time, I start the simulation with the companion just outside the envelope
from 124 R⊙, and injected energy below the surface with a half-Gaussian profile
peaking at the surface until the companion gets inside the envelope.

5. run_star_extras_agb_model4.f90 : Code for AGB model 4. I tried a simulation where I only take the ϕ component of the drag force on the
companion F2,ϕ. 

6. run_star_extras_agb_model5_input_3d_orbit.f90 : Model 5. In this case, I keep the orbit exactly the same as the 3D simulation, as well as the velocity. Goal was to see if we keep as many quantities same as the 3D simulation, what will be the envelope expansion. 





