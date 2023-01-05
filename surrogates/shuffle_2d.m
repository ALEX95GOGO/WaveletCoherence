function Surr=shuffle_2d(x,nSurr)

  %% shuffled surrogates of a data matrix x.

  %% INPUTS

    % x is data matrix
    % nSurr is the number of surrogates.

  %% OUTPUTS

    % Surr is the shuffled matrix

  [n_f,n_t] = size(x);
    for iSur = 1:nSurr
    % random vector
    % cell2mat(arrayfun(@(dim) randperm(dim), size(x), 'UniformOutput', false));
    rand_vec_f = randperm(n_f);
    rand_vec_t = randperm(n_t);
    % Shuffle the original x with respect to s1
    Surr(:,:,iSur) = x(rand_vec_f, rand_vec_t);
  end
end
