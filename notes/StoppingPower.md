## Fitting
gnuplot> a_c=.25;p_c=-.85;b_c=1.1;p2_c=.075                               
gnuplot> rep
gnuplot> fit c(x) 'Stopping power.YAG.jpg.collisional.dat' via a_c,p_c,b_c
iter      chisq       delta/lim  lambda   a_c           p_c           b_c          
0 6.0319896528e-01   0.00e+00  8.09e+00    2.500000e-01  -8.500000e-01   1.100000e+00
1 3.1040310875e-01  -9.43e+04  8.09e-01    2.515686e-01  -8.548533e-01   1.117750e+00
2 2.7034566683e-01  -1.48e+04  8.09e-02    2.583725e-01  -8.481571e-01   1.149831e+00
3 2.7027590192e-01  -2.58e+01  8.09e-03    2.591335e-01  -8.475576e-01   1.149893e+00
4 2.7027571627e-01  -6.87e-02  8.09e-04    2.590958e-01  -8.475904e-01   1.149925e+00
iter      chisq       delta/lim  lambda   a_c           p_c           b_c          

After 4 iterations the fit converged.
final sum of squares of residuals : 0.270276
rel. change during last iteration : -6.86857e-07

degrees of freedom    (FIT_NDF)                        : 19
rms of residuals      (FIT_STDFIT) = sqrt(WSSR/ndf)    : 0.119269
variance of residuals (reduced chisquare) = WSSR/ndf   : 0.014225

Final set of parameters            Asymptotic Standard Error
=======================            ==========================
a_c             = 0.259096         +/- 0.01642      (6.336%)
p_c             = -0.84759         +/- 0.01426      (1.682%)
b_c             = 1.14993          +/- 0.02799      (2.434%)

correlation matrix of the fit parameters:
a_c    p_c    b_c    
a_c             1.000 
p_c             0.994  1.000 
b_c            -0.510 -0.482  1.000 
gnuplot> rep
gnuplot> c(x)=a_c*(x**p_c)+b_c*(x**p2_c)      
