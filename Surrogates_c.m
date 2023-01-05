classdef Surrogates_c
    % class for Surrogate
    
    properties
        nSurr = 1,...  % number of surrogates
        maxiter = 1000,...  % maximum number of iterations allowed
        accerror = .001,...  % acceptable error
        error_change = 100,...  % acceptable relative error
        gradient = [0.8*pi pi],...  % amount where warped phase drops to determine parts of signale to shuffle
        p = 1,...  % order of the autoregressive part
        q = 0,...  % order of the moving average part
        ts = 0.1;  % sample time [s]
    end
    
    methods
    end
end

