


The program to run stuff in parallel is "distance_parallel" .  The usage is:

distance_parallel inputfile material_number epsilon volume_out pts_out num_processes

"volume_out" and "points_out" are the names it will use for output files.  points_out will be a text file with a list of points.  These files will be produced with different suffixes by the different processes, and the pts files must be cat'd together to make one big file.  Notice that "epsilon" is really epsilon, it will be squared before it is used.  So if you want epsilon sqrd to be 1.0e-12, then you make the call with epsilon set to 1.0e-6.

distance_parallel expects to find "distance3" in the same directory.  We can change this if we need to, but for now try to use a soft link to the executables in my directory.  

