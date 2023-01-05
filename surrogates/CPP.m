function Surr=CPP(ts_x,varargin)


%%  Cyclic Phase Permutation (CPP)

%% INPUT

  % ts_x time series matrix
  % nSurr (the number of surrogates you wish to generate)

%% OUTPUT

  % Surr the IAAFT surrogate time series

%% REFERENCES

  % adjusted CPP, original from 

Args = struct('nSurr', 1,...  % number of surrogates
              'gradient', [0.8*pi pi],...  % amount where warped phase drops to determine parts of signale to shuffle
                             ...  % larger gradient results in more conservative signi tests
              'disp_cycle_bins', false);  % plots cycles where phase drops

Args=parseargs_special(varargin,Args);
n_x = length(ts_x);

% Extract the instantaneous phase using the Hilbert transform
%   yh=hilbert(ts_x);
%   amplitude=abs(yh);
%   phase=(angle(yh));
%   phi = wrapTo2Pi(phase);
% Wrap angle in radians to [0 2*pi]
phi = wrapTo2Pi(ts_x);

% individual cycles by selecting points of discontinuity,
% here points where wrapped phase suddenly drops with pi
phase_diff = abs(phi(2:end)-phi(1:end-1));  % adjusted compaired to original algorithm (abs), because difference in both directions are relevant
small_diff = find(Args.gradient(1) >= phase_diff & phase_diff <= Args.gradient(2));

if ~any(small_diff)
    Args.gradient = Args.gradient*2;
    small_diff = find(Args.gradient(1) >= phase_diff &phase_diff <= Args.gradient(2));
end

%discontinuity_points = cell(length(small_diff)-1,1);

if Args.disp_cycle_bins
  figure;
  plot(phase_diff); hold on;
end

% should be optimized
for j = 1:length(small_diff)-1
    x_shuff = phi(small_diff(j)+1:small_diff(j+1));
    discontinuity_points{j} = x_shuff;
end

if any(small_diff)
    surr_start = phi(1:small_diff(1));
    surr_end = phi(small_diff(j+1)+1:end);
    %Surr = zeros(n_x, Args.nSurr);

    % Reconstruct a new phase series with these beginning and end points,
    % and randomly permutediscontinuity_points cycles
    for iSurr = 1:Args.nSurr
        Surr(:, iSurr) = unwrap(vertcat(surr_start,discontinuity_points{randperm(length(small_diff)-1)},surr_end));
    end
else
    warning('no phase shifts detected');
    Surr = ts_x;
end
