function varargout = surrogates(time_series,varargin)
  %% Surrogates
  %% REQUIREMENTS

    % Iterative Amplitude Adjusted Fourier Transform
    % Iterative Amplitude Adjusted Wavelet Transform
    % Cyclic Phase Permutation
    % Autoregressive Moving Average Process

  %% INPUTS

    % time_series: (length(time_series) by 1 data vector)
    % surr_name: type of surrogate algorithm: IAAFT, IAAWT, ARMA
    % nSurr: positive integer representing the number of surrogates
    % maxiter: positive integer representing is the  maximum number of iterations allowed of IAAFT
    % accerror: the acceptable error of IAAWT
    % error_change: relative error stop criterium of IAAWT
    % gradient: vector of two doubles, determindes where warped phase drops to determine parts of signale to shuffle
    % p: positive integer representing the order of the autoregressive part of the ARMA process
    % q: positive integer representing the order of the moving average part of the ARMA process
    % ts: is a double representing the sample time [s]

  %% OUTPUTS

    % Surr matrix of surrogates

Args=struct(    'surr_name', 'IAAFT', ...  % type of surrogate algorithm    
                'preprocessing', true);  % ensure
            
Args=parseargs_special(varargin,Args);
s = Surrogates_c; % import the class Surrogates_c

if Args.preprocessing
   input_ts = time_series;
  [time_series, ts_start, ts_end] = SurrPreprocessing(time_series, 'ts', s.ts);
end

n_x= length(time_series);

switch Args.surr_name
    case 'IAAFT'
          Surr = IAAFT(time_series,'nSurr', s.nSurr, ...
                        'maxiter', s.maxiter);
    case 'IAAWT'
          Surr= IAAWT(time_series,'nSurr', s.nSurr, ...
                        'accerror', s.accerror, ...
                        'error_change', s.error_change);
    case 'CPP'
          Surr= CPP(time_series,'nSurr',s.nSurr,...
                        'gradient', s.gradient);
    case 'ARMA'
          time_series=time_series(:); % Make time_series a column vector
          mean_ts = mean(time_series); % Mean of time_series
          std_ts = std(time_series);  % Standard deviation of time_series
          mean_ts = repmat(mean_ts,n_x,s.nSurr);
          std_ts = repmat(std_ts,n_x,s.nSurr);

          [alfa,beta,~,sigmaSquared] = ARMA2SR(time_series,s.p,s.q);
          err = sqrt(sigmaSquared)*randn(n_x,s.nSurr);
          beta = [1 beta'];
          alfa = [1 -alfa'];
          Surr = filter(beta,alfa,err);
          % Make surrogates have the same mean and variance of original series
          Surr = ProcessMatrix(Surr,0,1); % Normalize matrix of surrogates
          Surr = Surr.*std_ts+mean_ts; % Make surrogates have the same mean and
                                            % standard deviation as time_series
    otherwise
      warning('ERROR: not defined surrogate algorithm')
end
%% rebuild time series
if Args.preprocessing
  query = zeros(length(input_ts), size(Surr,2));
  query(ts_start:ts_end, :) = Surr;
  Surr = query;

end

varargout={Surr};
varargout=varargout(1:nargout);
end
