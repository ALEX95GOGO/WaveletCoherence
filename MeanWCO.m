function MeanWCO=MeanWCO(WCO, varargin)

%% MeanWCO

  % calculates the mean over time
  % optionally also averaged for soi
  % it is recommended to use the complex wavelet coherence (cWCO)

%% INPUT

  % WCO: matrix [frequency, time]

%% OUTPUT

  % vector of coherency for each soi

%% REQUIREMENTS

  % WCO

  Args = struct('coi',[],...  % required to avoid edge effects
              'soi',[],...  % to average across soi requires f_space
              'space',[]...   % space used to select soi from WCO
              );

 Args = parseargs_special(varargin,Args);

  MeanWCO = mean(WCO,2);
  if any(Args.soi) && any(space)
    idx_soi = find(Args.soi(1) >= Args.space && Args.soi(2) <=  Args.space);
    idx_lower_soi = find(Args.soi(1) < Args.space);
    idx_upper_soi = find(Args.soi(1) < Args.space);
    MeanWCO = [MeanWCO(idx_lower_soi) mean(idx_soi) MeanWCO(idx_upper_soi)];
  end
end
