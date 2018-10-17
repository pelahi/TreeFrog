# Running TreeFrog on AHF catalogue

First TreeFrog needs to be compiled, please follow the instructions on the Compiling part of the TreeFrog README. Once that is done the mpi flag needs to be switched off and the flag that tells TreeFrog the ID's do not correspond to an index (as with AHF halo ID's) needs to be switched on. These flags are found in CMakeCache.txt in the build/ directory.

The mpi flag is:

	TF_MPI:BOOL=ON

This needs to be set to OFF and the halo ID flag is:

	TF_HALOIDNOTINDEX:BOOL=OFF

This needs to be set to ON. Once that it done cmake needs to be run again before TreeFrog can be re-compiled:

	cmake ..  
	make all
Now the sample configuration file this is the treefog.cfg located in the base TreeFrog directory. There are three options that will need to be updated, the first one is:
	
	Input_tree_format = 3

The next one is:

	Max_ID_Value

This will be need to be set to the total amount of all particles in the simulation. The final option is:

	Nsteps_search_new_links

This enables TreeFrog's multi-snapshot linking, where this set the number of snapshot to do the linking over. TreeFrog can now be run from the build/ directory by:

	./treefrog -C ../treefrog.cfg -i <filelist> -s <numsnapshots> -o <baseoutputfilename>

Where \<filelist> is a file containing the _particles files for each snapshot from the AHF output.


