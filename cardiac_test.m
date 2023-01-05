function varargout = cardiac_test(x, varargin)
  %% test for signal fragments: cardiac Tests
    % uses power spectral density to evaluate if the heart rate is detectable in x
  %% REQUIREMENTS

    % hamming window

  %% INPUTS

    % x time series matrix
    % FOI frequency range of interest
    % ts sampling interval, in seconds
    % disp_pdf boolean option visualization of power spectral density function

  %% OUTPUTS

    % result is a table that contains means and the state where 1 represents if the cardiac oscillation is present
    % fspec is a struct including mean, sd, state for the frequency bin_names
    %       low frequency (lf), respiration, heart rate, high frequency (hf),
    %       foi (frequency range of interest), control (for comparison of state)

Args=struct('FOI', [0.1 0.5],...  % frequency range of interest
            'ts', 0.1, ...  % sampling interval, in seconds
            'disp_pdf', false);

Args=parseargs_special(varargin,Args);

%% Definition of frequency bins
fspec.lf.bin = [0 0.1];
fspec.respiration.bin  =[0.19 0.21];
fspec.heartrate.bin  = [1 1.6];
fspec.hf.bin  = [1.7 3];
fspec.foi.bin  = Args.FOI;
fspec.control.bin  = [0.19 1];

% number of time series
n = size(x,2); 
% number of points for discrete Fourier transform
nfft = max(256,2^(nextpow2(length(x))));

% pxx is the  power spectral density (PSD) estimate, returned as a real-valued,
% nonnegative column vector or matrix
% each column of pxx is the PSD estimate of the corresponding column of x
% f is a frequency vector in cycles per unit time == Hz
[pxx, f] = periodogram(x,hamming(length(x)),nfft,1/Args.ts);

if Args.disp_pdf
    figure;
    periodogram(x,hamming(length(x)),nfft,1/Args.ts); % or rect
end


bin_names=fieldnames(fspec);
result = table;
for bin=1:numel(bin_names)
      fspec.(bin_names{bin}).f_id = find(f > fspec.(bin_names{bin}).bin(1) & ...
                                f < fspec.(bin_names{bin}).bin(2));
      fspec.(bin_names{bin}).f = pxx(fspec.(bin_names{bin}).f_id,:);
      fspec.(bin_names{bin}).('mean')= mean(fspec.(bin_names{bin}).f,1)';
      fspec.(bin_names{bin}).('sd')= std(fspec.(bin_names{bin}).f,1)';
      query = table(fspec.(bin_names{bin}).('mean'), 'VariableNames',{bin_names{bin}});
      result = horzcat(result,query);
end

result.state=zeros(n,1);
if (sum(result.heartrate > result.control & result.heartrate > result.hf) > 0)
    result.state(result.heartrate > result.control & result.heartrate > result.hf) = 1;
end


varargout={result, fspec};
varargout=varargout(1:nargout);
end
