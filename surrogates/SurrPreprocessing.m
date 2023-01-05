function varargout=SurrPreprocessing(time_series,varargin)

  %% REFERENCES

  %  Wiener-Khinchin theorem
  % [1] Gemma Lancaster, Dmytro Iatsenko, Aleksandra Pidde,Valentina Ticcinelli, Aneta Stefanovska
  %     Surrogate data for hypothesis testing of physical systems Physics Reports (2018)
  % [3] C. J. Stam, J. P. M. Pijn, W. S. Pritchard,
  %    Reliable detection of nonlinearity in experimental time series with strong periodic components
  %   Physics Reports (1998) 361–380.

  Args=struct(    'amount_start', 0, ...  % proportion of time_series consider at the beginning
                  'amount_end', 0,...   % proportion of time_series consider at the end
                  'process_points', 0,...  % 10 recommanded [3] smoothing attention should be greater for noisy signals
                  'ts', 0.1);  % sample time [s]

  Args=parseargs_special(varargin,Args);

  time_series = time_series-mean(time_series);  % center time series
  n_t =length(time_series);


  if Args.process_points <= 0
    Args.process_points = 1/Args.ts;  % should be close to the sampling frequency
  end
  if Args.amount_start <= 0
    Args.amount_start = round(n_t/100);
  end
  if Args.amount_end <= 0
    Args.amount_end = round(n_t/10);
  end


%% Correct any mismatch between start and end points and corresponding first derivatives
  % finite mismatch between the beginning and end of the data, spurious results will be obtained
  % reused part of original time series
  prop_start = time_series(1:Args.amount_start);
  prop_end = time_series(end-Args.amount_end:end);

  % Truncate to match start and end points and first derivatives
  derivatives=zeros(length(prop_start)-Args.process_points,length(prop_end)-Args.process_points);

  for j=1:length(prop_start)-Args.process_points
      for k=1:length(prop_end)-Args.process_points
          derivatives(j,k)=sum(abs(prop_start(j:j+Args.process_points)-prop_end(k:k+Args.process_points)));
      end
  end
  derivatives = abs(derivatives);

  % Minimum mismatch
  [min_row,min_col]= find(derivatives == minmin(derivatives));

  ts_start = min_row;
  ts_end = min_col+length(time_series(1:end- Args.amount_end));
  Time_Series = time_series(ts_start:ts_end); % New truncated time series

  varargout={Time_Series, ts_start, ts_end};
  varargout=varargout(1:nargout);
end
    