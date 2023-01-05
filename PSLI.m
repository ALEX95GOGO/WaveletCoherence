function varargout = PSLI(wave_x,wave_y, periods, varargin)
%% Phase Slope Index (PSI)

 % Phase slope index (PSI) is implemented to estimate the time-lagged and linear interactions.
 % This method is done by investigating the complex-value coherency, and quantify the consistency of the 
 % direction of change in the phase difference(slope of the phase of cross-spectra) across frequency of interest.
 % A positive PSI will indicate the channel to be a driver and a negative PSI will indicate the channel to be a recipient.

%% INPUT

  % wave_x, wave_y matrix representing the spectra after wavelet transform of the signal 
  % periods is a vector , representing the period length for choosing the frequency band of interest
  % an option to normalization, 0 as no normalization and 1 as with normalization

%% OUTPUT

  % the phase slope index values. For a psi(i,j)>0, 
  %     it indicates the information flow from channel i to channel j.
%% REFERENCE

    % DOI: 10.1103/PhysRevLett.100.234101
    %      10.1016/j.neuroimage.2012.09.036 
    %      10.3389/fnsys.2015.00175


Args=struct('FOI', [0.12 0.8], ... %frequencey range of interest
            'normalization' , []); %0 for no normalization, 1 for normalization
                                 
Args=parseargs_special(varargin,Args);
if isempty(Args.normalization)
    Args.normalization = 0; %default as no normalization
end
%example: psi = PSLI(wave_x,wave_y, periods, 'normalization',1)
%% extract the frequency of interest
nchan = 2;  % number of channels
[nep,maxfreqbin, ~] = size(wave_x); % segment length, maximum frequency
spectra = zeros(maxfreqbin,2); 
periodOI = 1./Args.FOI;  % period of interest
indexRange = find(periods<periodOI(1)&periods>periodOI(2));

%% calculate the cross spectra
cs = zeros(nchan,nchan,length(indexRange));
i = 1;
for j = indexRange(1) : indexRange(end)  %choose the frequency band from the FOI
    spectra(:,1) = wave_x(j,:,1);
    spectra(:,2) = wave_y(j,:,1);

      
         for f=indexRange(1:end)  % for all frequencies
                 cs(:,:,i)=cs(:,:,i)+conj(spectra(f,:)'*spectra(f,:)); 
                 i = i+1;
         end
end

%%  Phase slope index calculation
     
df = 1;
[~, ~, nf]=size(cs);
pp=cs;
for f=1:nf
    pp(:,:,f)=cs(:,:,f)./sqrt(diag(cs(:,:,f))*diag(cs(:,:,f))'); % the complex coherence
end
psi=sum(imag(conj(pp(:,:,1:end-df)).*pp(:,:,1+df:end)),3);  % calculate the phase slope index
psloc = imag(conj(pp(:,:,1:end-df)).*pp(:,:,1+df:end));
stdpsi=squeeze(std(psloc,0,3))*sqrt(nep);
if (Args.normalization == 1)
    psi = psi./(eps + stdpsi);
end
%% output
varargout={psi,stdpsi};
varargout=varargout(1:nargout);
end