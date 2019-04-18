###vertical_crack
Details

    Horizon = 8 mm,  h = 2 mm, vertical crack of length 20 mm
    Displacement boundary condition
    Time = 100 ms, Delta t = 0.004 ms
    Number of threads = 4
Old code

    central_difference
        Memory: 12 MB
        Time: 69.3742 sec
    velocity_verlet
	    Memory: 12 MB
	    Time: 72.0021 sec	
New code

    central_difference
        Memory: 6.3 MB
        Time: 64.9109 sec
    velocity_verlet
        Memory: 6.3 MB
        Time: 68.0985 sec

###diagonal_pulling
Details

    Horizon = 4 mm,  h = 1 mm, angled crack of length 20 mm
    Force boundary condition
    Time = 140 ms, Delta t = 0.004 ms
    Number of threads = 4
Old code

    central_difference
        Memory: 21.7 - 24.5 MB
        Time: 412.3 sec
	velocity_verlet		
	    Memory: 21.7 - 24.5 MB
	    Time: 413.074 sec	
New code

    central_difference
        Memory: 8.7 - 14 MB
        Time: 357.874 sec
    velocity_verlet
        Memory: 8.7 - 14 MB
        Time: 362.198 sec
        
###diagonal_pulling_refine
Details

    Horizon = 2 mm,  h = 0.5 mm, angled crack of length 20 mm
    Force boundary condition
    Time = 140 ms, Delta t = 0.004 ms
    Number of threads = 4
Old code

	velocity_verlet		
	    Memory: 70.8 - 75.3 MB
	    Time: 1691.46 sec
New code

    velocity_verlet
        Memory: 27.3 - 45 MB
        Time: 1379.21 sec
			
##Summary
- New code takes smaller time consistently
- New code consumes less memory

##Other observation
- Each call to output function in new code increases the memory by some amount 
