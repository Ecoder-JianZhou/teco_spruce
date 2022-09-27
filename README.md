
# TECO SPRUCE 1.0


This repository includes fortran code to do data assimilation and forecasting responses of a northern peatland ecosystem to atmospheric CO2 fertilization and a gradient of experimental warming over next decade in the SPRUCE experiment (Jiang et al. 2018 JGR-Biogeosciences). 

TECO SPRUCE 1.0 can be run in two different environments either online or offline in desktop. To run the code online, please visit a web-based application Ecological Platform for Assimilation of Data (EcoPAD) at http://ecolab.cybercommons.org/ecopad_portal, which automates ecological forecasting at weekly timescale in SPRUCE. In brief, EcoPAD automates data transfer from SPRUCE ftp site, constrain parameters via data assimilation, and then forecasting ecosystem responses. The detail framework of EcoPAD and about installing in local machine is available at https://github.com/ou-ecolab/ecopad_documentation.  

Here, we walk through step by step to run the code offline in a desktop.

1.	Input data available for TECO SPRUCE 1.0

SPRUCE_pars.txt : initial/default parameter values for TECO model

SPRUCE_forcing.txt: forcing drivers include air temperature (Tair), soil temperature (Tsoil), relative humidity (RH), vapor pressure deficit (VPD), precipitation (Rain), wind speed (WS), and PAR.

SPRUCE_da_pars.txt: this file is only used for data assimilation. It contains three blocks of parameter values. The first block has values of “0” or “1”. “1” represents the parameter is chosen for data assimilation; “0” represents the parameter is not chosen for data assimilation, default values would then be used in simulations. The second and third blocks of parameter are minimum and maximum values, respectively.

SPRUCE_obs.txt: The pretreatment datasets from 2011to 2014 used for data assimilation. “days” means the number of days since 2011-01-01; “-9999” means the data is not available. See the paper for details on how the datasets were complied. 

2.	An empty “output” folder, which record simulation outputs, should be created. Here we uploaded three sample output files:

Paraest.txt: The constrained parameters from data assimilation, which can be used as input for simulation or forecasting by overriding the default values. The first row in the .txt file is a number. For example, here, “18” means 18 of the parameters were constrained by data assimilation. The second row is a series of numbers, which means the location of parameters in the default parameter pool (SPRUCE_par.txt). Then, the accepted parameter values were recorded as sequence.

Simu_dailyflux.txt: This is daily output from simulation of 2011-2014. The actual code to generate this file is located at lines 1256-1262 in the fortran code TECO_SPRUCE.f90. Users can modify the file to output any variables they are interested.

SPRUCE_yearly.txt: This is yearly output from the simulation. The actual code to generate this file is located at lines 1243-1251 in the fortran code TECO_SPRUCE.f90. Users can modify the file to output any variables they are interested. 
  
3.	Compile the fortran code TECO_SPRUCE.f90 (we assume users have basic fortran programing experience, and have complier installed)

First, compile the fortran code: gfortran TECO_SPRUCE.f90 –o TECO.out

Second, run the complied file using different argument options:
	
	 Option “0” model simulation (MCMC = 0), run the script:

 		 ./TECO.out input/SPRUCE_pars.txt input/SPRUCE_forcing.txt input/SPRUCE_obs.txt output/ 0 input/SPRUCE_da_pars.txt

	 Option “1” data simulation (MCMC = 1), run the script:

  		./TECO.out input/SPRUCE_pars.txt input/SPRUCE_forcing.txt input/SPRUCE_obs.txt output/ 1 input/SPRUCE_da_pars.txt


	 Option “2” forecasting (MCMC = 2), run the script:

		./TECO.out input/SPRUCE_pars.txt input/SPRUCE_forcing.txt input/SPRUCE_obs.txt output/ 2 input/SPRUCE_da_pars.txt input/Weathergenerate 2024 365 2.25 380.0

The detail argument format and description is uploaded to argument_format.txt file.

