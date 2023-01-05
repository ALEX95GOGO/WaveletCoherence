load('HS_Mit_M_N856_MES_Probe1.mat')
ts_x1 = filtert.hbo(:,8);
load('HS_Mit_M_N856_MES_Probe2.mat')
ts_x2 = filtert.hbo(:,8);

%% WCOExperiment Testing
parameters.wavelet_type = 'gmw';
parameters.be = 3;
parameters.ga = 3;
parameters.wave_normalization = 'energy';
parameters.ws_size = 8;
parameters.wt_size = 8;
parameters.window_time = 'hann';
parameters.surr_name = 'ARMA';
parameters.coherence = 'WCO';
parameters.signilevel = 0.01;
parameters.smoothing = 'Aguiar';
parameters.preprocessing = 0;
WCOExperiment(ts_x1, ts_x2, parameters);
