classdef LocalInjectedSynchronySimulation
    % Class for Local Injected Synchrony Simulation (LISS)
    
    properties
         FOI = [0.078 0.495];  % frequencey range of interest
         ts = 1/10;  % ts sampling time [s]
         T = [120, 200, 300];  % periodlength of short RECT (see HDR duration)
         tr = [1000, 1500, 2000];  % translatation of short RECT
         shift_time = 0;  % integer, shifts the simulated synchrony in time
         shift_freq = 0;  % double, shift in frequency domain
         snr = 1;   % singnal to noise ratio, integer, for <=1 no noise will be added
                     ...  % > 1 gaussian white noise will be added to the manipulated signal Y
    end
    
    methods
    end
end

