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

4. run_star_extras_agb_model3.f90 : Code for AGB simulation Model 3. In the first part of the simulation, up to 95 days, I neglect the motion of the
core with respect to the envelope. Consequently, by ignoring the motion of the core, I only calculate the ϕ compo-
nent of the drag force on the companion(F2,ϕ). From 95 days onwards, I again take into account the motion of the core.
Accordingly, I take the drag force to be Fϕ, calculated with the same method
as in model 1. To solve the time mismatch issue, instead of reading the drag force data from
a later time, I start the simulation with the companion just outside the envelope
from 124 R⊙, and injected energy below the surface with a half-Gaussian profile
peaking at the surface until the companion gets inside the envelope.

5. run_star_extras_agb_model4.f90 : Code for AGB model 4. I tried a simulation where I only take the ϕ component of the drag force on the
companion F2,ϕ. 

6. run_star_extras_agb_model5_input_3d_orbit.f90 : Model 5. In this case, I keep the orbit exactly the same as the 3D simulation, as well as the velocity. Goal was to see if we keep as many quantities same as the 3D simulation, what will be the envelope expansion. 





