function Surr=shuffle_1d(x,nSurr)

  %% shuffled surrogates of a data matrix x.

  %% INPUTS

    % x is data vector
    % nSurr is the number of surrogates.

  %% OUTPUTS

    % Surr is the shuffled vector.

   [n_t, n_f] = size(x);
   if n_f > 1
       warning('shuffle_d1:false dimensions');
   end
   for iSur = 1:nSurr
    % random vector
    rand_vec_t = randperm(n_t);
    % Shuffle the original x with respect to s1
    Surr(:,iSur) = x(rand_vec_t, 1);
  end
end
