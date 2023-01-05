function [Surr,Iter]=IAAFT(x,varargin)

%% Iterative Amplitude Adjusted Fourier Transform (IAAFT)


%% INPUTS

  % x (length(x) by 1 data vector)
  % nSurr (the number of surrogates you wish to generate)
  % maxiter is the  maximum number of iterations allowed. maxiter must be scalar

  %% OUTPUTS

  % Surr the IAAFT surrogate time series
  % Iter is the number of iterations needed of the i-th surrogate series

%% REFERENCES

  % Adapted from Dr D. Kugiumtzis & Alexandros Leontitsis: Chaotic Systems Toolbox
  % Schreiber, Thomas, Schmitz, Andreas "Improved surrogate data for nonlinearity tests" Phys. Rev. Lett. 77, 635 (1996)

Args = struct('nSurr', 1,...  % number of surrogates
            'maxiter', 1000);  % maximum number of iterations allowed

Args = parseargs_special(varargin,Args);

% The magnitudes of x
amp = abs(fft(x));

% Shuffle x
Surr = shuffle_1d(x,Args.nSurr);

% Sort x
[x,r] = sort(x);


for j=1:Args.nSurr

    % Calculate the phases of the shuffled series
    phase = angle(fft(Surr(:,j)));

    % Initialize the loop
    k = 1;
    indold = r;
    converge = 0;
    while k<=Args.maxiter & converge == 0
        % Make phase-randomized surrogates ...
        Surr(:,j) = amp.*exp(phase.*i);
        Surr(:,j) = real(ifft(Surr(:,j)));
        % ... and give them the distribution of x
        [Surr(:,j),T] = sort(Surr(:,j));
        [Surr(:,j),indnew] = sort(T);
        Surr(:,j) = x(indnew);
        % Check the convergence
        if indnew==indold
            converge = 1;
        else
            indold = indnew;
            k = k+1;
        end
        % Loop again if needed, calculating the phases once more
        phase = angle(fft(Surr(:,j)));
    end

    % Get the iterations of each surrogate series
    Iter(j) = k;

end
