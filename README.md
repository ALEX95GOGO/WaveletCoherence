# WaveletCoherence

This toolbox provides different algorithms to analyze the wavelet spectra of time series signals and the wavelet coherence for wavelet spectra.

Functions:

**CWT.m**     
**WCO.m**     
**surrogates.m**   
**cardiac_test.m**   

**CWT**

The *CWT.m* provides the continuous wavelet transformation([1],2.4) of a single time series or a matrix of multiple time series. The function supports two different types of wavetlets: the generalized morse wavelet(gmw)([1],2.9.2) and the morlet wavelet([1],2.9.1). The implementation relies mainly on the Matlab cwt.m fuction as well as the the JWavelet module of [JLab toolbox](http://www.jmlilly.net/doc/jWavelet.html). 

  * supports different wavelets (gmw, morlet wavelet)    
  * wavelet parametization   
  * wavelet normalization options   
  * supports frequency space in scale or period length  
  * different calculation of the frequency resolution  
  * supports real as well as imaginary time series  
  * input of single as well as multiple time series  
  * calculate the cone of influence (coi)  
  * plot the wavelet spectrum as well as the time series, coi and space of interest (area between two white lines)   
  * plot the mother wavelet   
  * returns the heisenberg area, radius(standard deviation) in time and frequency space of the gmw   
*Features*:  

**WCO**

The *WCO.m* provides the wavelet coherence([2],2.4) of two wavelet spectra. The shannon entropy([3],pg.956) and phase locking value([4],pg.7) of the wavelet coherence can also be calculated by this function. The implementation relies mainly on the [wavelet coherence toolbox](http://grinsted.github.io/wavelet-coherence/) by Aslak Grinsted, the [ASToolbox](https://sites.google.com/site/aguiarconraria/joanasoares-wavelets/the-astoolbox), and the WCOH toolbox by Cohen[5].

*Features*:

  * support first order wavelets as well as higher-ordered multiwavelet   
  * option for output real wavelet coherence or complex wavelet coherence
  * smoothing parametization    
  * support different smoothing algirithm   
  * plot the wavelet coherence spectrum as well as the time series, coi and space of interest(area between the white lines)   
  * returns the wavelet coherence, phase locking value and shannon entropy  

**surrogates**

The *surrogates.m* provides surrogate time series through iterative amplitude adjusted wavelet transform.

*Features*:  
  * supports different surrogate algorithm   
  * iterations parametization   
  * supports multiple numbers of surrogates   

**Cardiac_test**

the *cardiac_test.m* evluate if the heart rate is detectable in a time series using the power spectral density. If the state = 1 in result, the cardiac oscillation is present.   

*Features*:  
  * supports time series matrix(multiple time series)  
  * plot the power spectral density function  
  * shows mean and sd for different frequency bins (low frequency (lf),respiration, heart rate, high frequency(hf), foi(frequency of interest) and control.

*Reference*:   
[1] Aguiar-Conraria, L. and Soares, M.J. (2010) "The Continuous Wavelet Transform: A Primer", NIPE Working paper   
[2] Cazelles, Bernard, et al. "Time-dependent spectral analysis of epidemiological time-series with wavelets." Journal of the Royal Society Interface 4.15(2007):625.   
[3] Cazelles, Bernard, and L. Stone. "Detection of imperfect population synchrony in an uncertain world." Journal of Animal Ecology 72.6(2010):953-968.   
[4] Bastos, A. M., and J. M. Schoffelen. "A Tutorial Review of Functional Connectivity Analysis Methods and Their Interpretational Pitfalls. " Frontiers in Systems Neuroscience 9.Pt 2(2016):175.    
[5] Cohen, Ed A. K. and Andrew T. Walden. "A Statistical Analysis of Morse Wavelet Coherence."IEEE Transactions on Signal Processing 58 (2010): 980-989.    