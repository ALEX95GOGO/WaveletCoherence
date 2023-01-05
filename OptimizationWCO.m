

%% REFERENCE
  % [2] https://de.mathworks.com/help/stats/bayesian-optimization-algorithm.html#bvaz8tr-1
  % Jasper Snoek, Hugo Larochelle and Ryan P. Adams (2012)
  % "PRACTICAL BAYESIAN OPTIMIZATION OF MACHINE LEARNING ALGORITHMS"


%% start data
load('test_cwt_equal_JoVe.mat');
ts_x1 = filtert_p1.hbo(:,1);
ts_x2 = filtert_p2.hbo(:,12);

%% simulation szenarios
% [ts_x1, ts_x2, labels_vec, rect_vec, labels_shifted_vec] = LISS(ts_x1,ts_x2, ...
%             'filter','bandstop',...
%             'shift_time', 0,...
%             'shift_freq', 0.1,...
%             'snr',1,...
%             'T', [100, 110, 120, 200, 200, 500],...  % periodlength of short RECT (see HDR duration)
%             'tr', [500, 800, 1100, 1500, 2000, 2800],...  % translatation of short RECT
%             'disp_info',true);
[ts_x1, ts_x2, labels_vec, rect_vec, labels_shifted_vec] = LISS(ts_x1,ts_x2, ...
            'filter','bandstop',...
            'shift_time', 10,...
            'snr',1,...
            'T', [ 110, 120, 120],...  % period length of short RECT (see HDR duration)
            'tr', [1200, 2000, 2800],...  % translatation of short RECT
            'disp_info',true);


%% WCOExperiment Testing
parameters.wavelet_type = 'gmw';
parameters.be = 20;
parameters.ga = 3;
parameters.wave_normalization = 'energy';
parameters.ws_size = 5;
parameters.wt_size = 5;
parameters.window_time = 'hamming';
parameters.surr_name = 'ARMA';
parameters.coherence = 'WCO';
parameters.signilevel = 0.05;
parameters.smoothing = 'Aguiar';
parameters.preprocessing = 0;
parameters.morse_space = 'ag_manuel';

%WCOExperiment(ts_x1, ts_x2, labels_vec, parameters);
%% Parameters
be = optimizableVariable('be',[3,20],'Type','integer');
ga = optimizableVariable('ga',[3,6],'Type','integer');
wave_normalization = optimizableVariable('wave_normalization',{'energy','bandpass'},'Type','categorical');
ws_size = optimizableVariable('ws_size',[5,10],'Type','integer');
wt_size = optimizableVariable('wt_size',[5,10],'Type','integer');
surr_name = optimizableVariable('surr_name',{'ARMA','IAAFT', 'CPP', 'IAAWT'},'Type','categorical');
coherence = optimizableVariable('coherence',{'WCO', 'WCO2'},'Type','categorical');

signilevel =  optimizableVariable('signilevel',{'0.05', '0.01'},'Type','categorical');
smoothing = optimizableVariable('smoothing',{'Aguiar', 'Matlab', 'Cui'},'Type','categorical');  % + 'Cohen'
wavelet_type = optimizableVariable('wavelet_type',{'gmw', 'amor'},'Type','categorical');
window_time = optimizableVariable('window_time',{'Bartlett', 'blackman', 'hamming', 'hann'},'Type','categorical');  % here we assume that the window is id for time and scale
kappa = optimizableVariable('kappa',[15,22],'Type','integer');  % smoothing parameter for coen
morse_space = optimizableVariable('morse_space',{'density', 'length', 'matfreq', 'ag_manuel', 'ag_auto'},'Type','categorical');


%preprocessing = optimizableVariable('preprocessing',[0, 1],'Type','integer');  % logical values have to be integers
rng default  % set seed
%c = cvpartition(length(ts_x1),'Kfold',5);  % index for partitioning the set of time series

fun = @(x)OptimizationWCOHelper(ts_x1, ts_x2, labels_vec, x);


results = bayesopt(fun,[be,ga, wave_normalization, ws_size, wt_size, window_time, surr_name, coherence, signilevel, smoothing, wavelet_type, kappa, morse_space],...
        'Verbose',0,...  %  Command-line display off
        'AcquisitionFunctionName','expected-improvement-plus',...  %  the expected amount of improvement in the objective function more in [2]
        'MaxObjectiveEvaluations',10,...
        'UseParallel',false...  % requires
        )
  
