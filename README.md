LApx (Los Alamos Polycrystal code) predicts the mechanical response of crystalline materials subjected to arbitrary quasi-static (i.e., low strain rates) scenarios. 
LApx uses fast Fourier transform based methods [1] to solve for the stress and strain states at each point within a complex microstructure. 
Uniquely, this capability predicts the relative and absolute contributions of all plastic deformation mechanisms activated during loading. 
These mechanisms include dislocation glide, dislocation climb [2], vacancy diffusion processes (i.e., diffusion along grain boundaries and within the bulk, 
effect of non-equilibrium point defects due to irradiation) [3], and porosity evolution due to damage.

The code LApx solves the mechanical equilibrium condition in a voxelized microstructure. 
Key to it is the use of a Fast Fourier Transform based algorithm [1] to rapidly quantify the stress and strain fields at each point. 
Further, the code is equipped with several advanced constitutive equations which relate the inelastic strain at each point to both stress, 
temperature and features of the microstructure (i.e., dislocation content, precipitate content, local porosity) [3,4]. 
All of these features of the constitutive equations have been the subject of either reports or publications available in the open literature.

Key-Contributors and Developers:
Andrea Rovinelli
Laurent Capolungo
Ricardo A Lebensohn
Carlos Tomè

Bibliography
[1] Lebensohn, Ricardo A., Anand K. Kanjarla, and Philip Eisenlohr. 
"An elasto-viscoplastic formulation based on fast Fourier transforms for the prediction of micromechanical fields in polycrystalline materials." 
International Journal of Plasticity 32 (2012): 59-69.

[2] Lebensohn, Ricardo A., Craig S. Hartley, Carlos N. Tomé, and Olivier Castelnau. 
"Modeling the mechanical response of polycrystals deforming by climb and glide." 
Philosophical Magazine 90, no. 5 (2010): 567-583.

[3] Wen, Wei, A. Kohnert, M. Arul Kumar, Laurent Capolungo, and Carlos N. Tomé. 
"Mechanism-based modeling of thermal and irradiation creep behavior: An application to ferritic/martensitic HT9 steel." 
International Journal of Plasticity 126 (2020): 102633.


[3] Wen, Wei, A. Kohnert, M. Arul Kumar, Laurent Capolungo, and Carlos N. Tomé. 
"Mechanism-based modeling of thermal and irradiation creep behavior: An application to ferritic/martensitic HT9 steel." 
International Journal of Plasticity 126 (2020): 102633.

[4] Wang, H., Laurent Capolungo, Bjorn Clausen, and C. N. Tomé. 
"A crystal plasticity model based on transition state theory." 
International Journal of Plasticity 93 (2017): 251-268.
