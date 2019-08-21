This Master Thesis project was done by 
Zoltán Márk Pintér

 (s172040)
Email: pinterzoltanmark@gmail.com
In collaboration with Danfoss A/S RAC

Advisor(s):
Hans Henrik Niemann,
Dimitrios Papageorgiou,
Lars Finn Sloth Larsen (external)



DTU Elektro

Technical University of Denmark
2800 Kgs. Lyngby

Denmark


elektro@elektro.dtu.dk 


The most relevant content are the following. The rest are support scripts for repetitive tasks, for plotting and saving, and unused codes of the special course.

------------------------------------------------------	CONTENT	-----------------------------------------------------------

Special_Course_on_Modelling_a_Refrigeration_System.pdf  	Documentation of the special course
Codes 								Folder containing codes
 - Linearization 	 					Folder including linearization related files
 -  - analysis.m 	 					Linear time invariant analysis (stability, obsv, ctrb)
 -  - KalmanFilter.m  						Kalman filter object with predictor
 -  - linearization.m  						Running the jacobian functions
 -  - main_lin.m 	 					Main file to run the codes in the proper order
 -  - modelcreation.m  						Creation the original complex model
 -  - modelcreation_simplified.m	 			Creation of the final simplified model
 -  - normalization.m 						Normalizing the model for analysis
 -  - partialdifferentials.m					Visualizing the partial differentials
 -  - propertyconversion.m 					Finding the LSQ transfer parameters
 -  - substitution.m  						Substituting values for LTI analysis, complex system
 -  - substitution_simplified.m 				Substituting values for LTI analysis, simplified system
 - MatlabPlant 							Folder including the scripts and objects of the special course
 -  - PIController.m 						PI controller object
 - Measurements 						Folder including data and codes related to the field experiment
 -  - findstruct.m 						Finding the structures of the smallest test error for every model order
 -  - Gasloop_20190606.csv 					Field experiment data
 -  - main_meas.m 						Running the estimatortest.m file for the field data
 -  - paramanal.m 						Analysis of found parameters 
 -  - readcsv.m 						Reading the field experiment data
 -  - residualwhitening.m					Main file for whitening the residual 
 -  - structbuilder.m 						Finding the best fit for all structures 
 - ModelicaPlant 						Folder including the main part of the thesis work
 -  - bestsimulation 						Folder for storing the most descriptive simulation data
 -  - matfiles 							Folder for rest of the matfiles 
 -  - object_inits 						Folder for the objects of this thesis and their initializations
 -  - - FaultDetector.m 					Fault Detector object for GLR, CUSUM and EM
 -  - - RecursiveLeastSquares.m 				RLS object
 -  - simulink 							Folder for storing the simulink files
 -  - - AmbientTester.slx 					Ambient temperature design file
 -  - - DistributedControllersWithFault.slx 			Original stress tests for detectors and estimators
 -  - - Estim_DisCont_Fault.slx 				Final simulation generator
 -  - - SupermarketPACK.fmu 					Mathematical simulation file transferred from Dymola			
 -  - c2dn.m 							Time discretization
 -  - constantsave.m 						Saving of constant system parameters
 -  - constrainer.m  						Function for applying steady-state constraint
 -  - crosscov.m 						File for finding the virtual air temperature time constant 
 -  - distributedControllers.m 					Discretization of continuous time controllers designed in Modelica
 -  - estimatortest.m  						Script to stress-test the estimator
 -  - estimator_simulink.m 					Script for evaluating results from Simulink simulations
 -  - g.m 							Nonlinear output function
 -  - GCstatemonitoring.m 					File for investigating the enthalpy or temperature states within the gas cooler
 -  - LTVsystemDescription.m 					Setpoint selector function
 -  -  - x2q.m 							Conversion from receiver liquid level to vapour quality