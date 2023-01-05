%% Additive White Gaussian Noise for Time Series
% adds white Gaussian noise to signal based on the signal to noise ratio
% calssical noise model used in Information theory
% additive: noise will get added to transmitted signal
% white: noise is equally present with the same power at all the frequencies
% gaussian: gaussian distribution with 0 mean and variance as the noise power
% for users without Communications System Toolbox

%% INPUT

  % ts_x: real valued time series vector
  % snr: signal to noise ratio, integer [1:] in dB, is used to determine the appropriate noise level

%% OUTPUT

  % ts_X: time series with added white gaussian noise

function ts_X = WhiteGaussianNoise(ts_x, snr)
    ts_n = length(ts_x);
    snr = 10^(snr/10);
    energy = sum(abs(ts_x).^2)/(ts_n);  % determine the energy
    nsd = energy/snr;  % noise spectral density for same power for all the frequencies
    %% calulate a random noise regarding the snr
    noiseSigma = sqrt(nsd);
    noise = (noiseSigma*randn(1,ts_n))';  % random vector using the standard deviation of the noise spectral density

    ts_X = ts_x + noise;  % add noise to input signal --> additive noise
end
