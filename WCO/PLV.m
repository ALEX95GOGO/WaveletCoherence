%% phase locking value
function [Phi,Entropy] = PLV (wco)
  Phi = atan (imag(wco)./real(wco));
  Entropy = wentropy(wco,'shannon');
end
