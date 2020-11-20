
########################################################
##                                                    ##
##   Instructions for training auto grid regression   ##
##   model for 2+1d nonconformal hydrodynamics        ##
##                                                    ##
########################################################


1) Generate training data (e.g. anisotropic hydro, 1000 model parameter samples, 1 smooth Trento event)


	sh generate_training_data_smooth.sh vah 1000 1		cpu_vah/scripts/auto_grid/generate_training_data

	args:	$1 = hydro mode
			$2 = number of training parameter samples (determines number of jobs)
			$3 = number of hydro events per job for this training session

	(wait until all the jobs have finished...)




2) Process training data and fit the regression model


	sh fit_regression_model.sh vah 1000 1				cpu_vah/scripts/auto_grid

	note: commands arguments are the same as step 1


3) Generate new model parameter samples and predict the fireball radius (40 new model parameter samples)


	sh predict_fireball_radius.sh vah 40 1				cpu_vah/scripts/auto_grid

	args:	$1 = hydro mode
			$2 = number of new parameter samples (for launch)
			$3 = number of hydro events per job during training (see step 1)

	note: $1 and $3 are the same as steps 1 and 2
	note: model_parameters and fireball_size_predictions are stored in cpu_vah/python/




4) Launch regression model (e.g. run 1 hydro event with model_parameters_6.dat)


	set auto_grid = 1 									cpu_vah/parameters/lattice.properties

	sh hydro.sh 1 6 									cpu_vah (max value of 2nd command argument = number of new samples generated)

	note: auto grid isn't used for default model parameters
	note: make sure you run the correct hydro mode (i.e. vah)




########################################################
##                                                    ##
##   Optional scripts after training the auto grid    ##
##                                                    ##
########################################################


1) Benchmark auto grid against large fixed grid (anisotropic hydro, 50 fluctuating events per parameter samples)

	sh benchmark_fixed_grid.sh vah 50 					cpu_vah/scripts/auto_grid

	args:	$1 = hydro mode
			$2 = number of hydro events per job for this benchmark test

	note: $1 is the same as step 1
	note: number of test model parmeters is fixed to 200

	(wait until all the jobs have finished...)


	sh benchmark_auto_grid.sh vah 50 1 2.5 0 			cpu_vah/scripts/auto_grid

	args:	$1 = hydro mode
			$2 = number of hydro events per job for this benchmark test
			$3 = number of hydro events per job during training (see step 1)
			$4 = margin parameter
			$5 = sigma factor

	note: $1 and $3 are the same as step 1

	(wait until all the jobs have finished...)


	Jupyter notebook: benchmark_runtimes.ipynb			cpu_vah/python




2) Visualize the training data

	Jupyter notebook: visualize_training_data.ipynb		cpu_vah/python









