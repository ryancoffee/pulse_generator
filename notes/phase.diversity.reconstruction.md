## Phase Diversity Reconstruction  
Inspired by talk from   
Serge Bielawski  

Han et al., IEEE Tansactions in Microwave Theory and Technology v53 p1404 y2005  
Transfer functions  
Limitation on temporal resolution from chirp versus window ~ sqrt(window * Fourier Limit)  
from Sun et al., Apllied Phys Lett 73 2233 1989  

But *this is no longe true* if we use say two chirps simultaneously, or so he says the two polarization states, simultaneously  
Then we are basically like measuring in quadrature, so there should be no zeros in the transfer functions  
(H1(X) * Y1 + H2(X) * Y2) / (H1(X)^2 + H2(X)^2)  
Something like this, we see that there should be no zeros in the denom.

## Neural Networks   
Basically the hypothesis is that the zeros from the transfer function are what kill a deconvolution  
Can we treat this as an inverse problem that we could solve with a trainied neural network  
If we use "Phase Diversity Reconstruction" then we can do it... but if we cannot, as in the PET case likely...  
Can we compensate for this by training to solve the inverse problem?  
