classdef ComplexWaveletTransform
    % class for COMPLEX WAVELET TRANSFORM (CWT)
    
    properties
        freq_limits = [0.0156, 5]; % low and high frequency boundaries used to calculate the frequency bins (morse space and morlet space)
        ts = 0.1; % sampling interval, in seconds
        n_v = 10; % number of voices per octave:also called frequency resolution, determines the number of samples range:[4-48] default 10
        K = 1; % the order of gmw
        be = 20; % beta for gmw
        ga = 3; % gamma for gmw
        energy_density = 4; % descibes the resolution in the morse space used for 'morse_space' == 'density'---> increases with overlap in frequency domain [REF]
        soi = [2.02 12.8]; % space of interest border (here for period length), unit must be the same as graph
                           ...  % should not violate nyquist criteria
    end

methods
    
    function [energy_freq, inst_freq, peak_freq]=gmw_peakf(obj)
        %  see Aguria
        
        if (obj.be ~= 0)
            peak_freq =  exp(1/obj.ga*(log(obj.be)-log(obj.ga)));
        else
            peak_freq = log(2)^(1/obj.ga);
        end
        
        energy_freq = exp( gammaln((2*obj.be+2)/obj.ga)-gammaln((2*obj.be+1)/obj.ga) )*2^(-1/obj.ga); %source from Aguiar, GMWMeasures
        inst_freq = exp(gammaln((obj.be+2)/obj.ga)-gammaln((obj.be+1)/obj.ga));
        
    end
end
end

