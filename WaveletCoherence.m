classdef WaveletCoherence
    % class for WAVELET COHERENCE(CWT)
    
    properties
        soi = [2.02 12.8];  % space of interest border (here for period length), unit must be the same as graph
        kappa = 20;  % smoothing parameter for cohen
        d = 1; % morlet wavelet shape diameter for cohen
        n_v = 10;  % number of octaves: this will do 10 sub-octaves per octave
        ts = 0.1; % sampling time [s]
        ws_size = 8; % to determine the size of the Hamming window window used for smoothing in the scale direction; size depends on dj and is  given by: ws_size/(2*dj) (with a minimum value of 5)
        wt_size = 8; % to determine the size of the Hamming window used for smoothing
                            % in the time direction; sizevaries with scale s and
                            % is given by:  wt_size*s/dt (with a minimum value of 5
    end
    
    methods
    end
    
end

