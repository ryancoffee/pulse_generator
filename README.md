# generate\_dataset
==============================
Generating dataset of broadband pulses of varying spectral phase and noise floor.

## Notes to self
-------------
For what its worth, this used up only 1/2 of the 32GB RAM on my laptop using 4 cores dual threaded for total 8 threads in parallel.
```bash
coffee@coffee-ThinkPad-X395:~/projects/slac/pulse_generator (main)$ top

top - 13:44:28 up  4:54,  1 user,  load average: 3.74, 1.09, 0.54
Tasks: 313 total,   3 running, 309 sleeping,   1 stopped,   0 zombie
%Cpu(s): 99.8 us,  0.2 sy,  0.0 ni,  0.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
MiB Mem :   5632.6 total,    138.8 free,   4725.5 used,    768.2 buff/cache
MiB Swap:   2048.0 total,   2019.2 free,     28.8 used.    620.6 avail Mem 

    PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                                                                                                     
  28823 coffee    20   0 3361504   2.7g   7904 R 792.4  49.2   4:24.75 generate_datase                                                                                                                             
   1640 coffee    20   0 1266276  79752  22800 S   2.3   1.4   7:55.35 Xorg                                                                                                                                        
   3146 coffee    20   0  826468  47380  27456 S   2.0   0.8   1:26.09 gnome-terminal-   
```

-------------
And the code timing was,  
```bash
coffee@coffee-ThinkPad-X395:~/projects/slac/pulse_generator (main)$ ./generate_dataset 
		======================================
		======= gen pulses started ========
		===== Tue Jun 13 13:43:53 2023
		===== on host coffee-ThinkPad-X395:::::
		===== using 8 threads
		======================================
filebase in main() is /home/coffee/pulsegendata/h5files/dataset
initializing pulse and plans
	It has taken 0 s for initializing pulse and fftw plans
waves.size()	65536	waves.front().size()	11000
	loop n = 	loop n = 	loop n = 00 in thread  in thread 6
	loop n = 0 in thread 4	loop n = 	loop n = 
	loop n = 70 in thread 0
0	loop n = 0 in thread 2
 in thread 
1
0 in thread 3
0 in thread 5




 ---- just left parallel region -----
 ---- and destroyed plan vectors ----
 ---------- runtimes are ------------
47	47	47	47	47	47	47	47	
 ------------------------------------
 ---------- Writing H5 file ---------
 --------- consider parallel --------
outfile = /home/coffee/pulsegendata/h5files/dataset-2023-6-13-h13-m44.h5
		======================================
		======== generate pulses stopped =======
		===== Tue Jun 13 13:44:58 2023
		===== in 65 s ====
		======================================

```


Project Organization
--------------------

    .
    ├── AUTHORS.md
    ├── LICENSE
    ├── README.md
    ├── makefile
    ├── config/
    ├── src/
    ├── include/
    ├── bin/
    ├── models/
    ├── notebooks/
    ├── sandbox/
    ├── docs/
    ├── figs/
    ├── references/
    ├── objects/
    └── reports/
