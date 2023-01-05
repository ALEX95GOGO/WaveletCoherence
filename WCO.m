function varargout=WCO(cwt_x,cwt_y,scales,periods,varargin)
%% Wavelet coherence
% provide different algorithms to calculate wavelet coherence for two wavelet spectrum
  %% REQUIREMENTS
  
    % Smoothing Operation for Cross Wavelet spectra
    % Wavelet Coherence
    % Phase Locking Value
    % Shannon entropy
    
  %% INPUTS
  
    % two wavelet spetrum
    % scales and period length
    % optional input for visulization (cone of influence, time series)
    % smoothing agirithm
    % coherence type
    % frequency of interest in herz
    % sampling time in seconds
    % smoothing parameter (for Cohen)
    % morlet wavelet shape diameter (for Cohen)
    % number of octave per wave
    % size of window in time or scale direction (for Aguiar)

  %% OUTPUTS
  
    % Wavelet coherence
    % Phase locking value
    % Shannon entropy
    
  %% References:
  
    %   [1] Aguiar-Conraria, L. and Soares, M.J. (2010)
    %       "The Continuous Wavelet Transform: A Primer", NIPE Working paper
    %   [2] Torrence, Christopher, and Peter J. Webster. "Interdecadal 
    %       changes in the ENSO�Cmonsoon system." Journal of Climate 12.8 
    %       (1999): 2679-2690.
    %   [3] Cohen, Ed AK, and Andrew T. Walden. "A statistical study of 
    %       temporally smoothed wavelet coherence." IEEE Transactions on 
    %       Signal Processing 58.6 (2010): 2964-2973.
    
    
%%
Args=struct('coi',[],...  % coen of influence recommended for plot
            'time_series',[],...  % time series recommended for plot
            'smoothing','Aguiar',...  % default as cohen, can choose windowing Torence or multiwavelets
            'coherence','cWCO',...  % default as cWCO, can choose WCO or WCO2
            'disp_coherence', false, ... %option to plot the coherence
            'y_axis','period',...  % can choose Scale and Period
            'window_time','hann',... % smoothing in time hamming, hann, blackman, Bartlett, otherwise no applied smoothing
            'window_scale','hann');  % smoothing in time hamming, hann, blackman, Bartlett, otherwise no applied smoothing


Args=parseargs_special(varargin,Args);
s = WaveletCoherence;
% required for window evaluation (preimplemented in Matlab)
supported_windows = {'Bartlett', 'blackman', 'hamming', 'hann'};

%% exception handling for higher order input
if (size(size(cwt_x),2) > 3 || size(size(cwt_y),2) > 3 )
    warning('WCO: input reaches critical dimensions!');
end

order = size(cwt_x,3);
if (order > 1 &&  ~strcmp(Args.smoothing, 'multiwavelets'))
  warning('WCO: smoothing changed to multiwavelets')
  Args.smoothing = 'multiwavelets';
end
switch Args.smoothing
%% Cohen's method
 %only for morlet wavelet
    case 'Cohen'
      [sWxx] = temporal_smooth(cwt_x,cwt_x,scales,s.kappa,s.ts,s.d);
      [sWyy] = temporal_smooth(cwt_y,cwt_y,scales,s.kappa,s.ts,s.d);
      [sWxy] = temporal_smooth(cwt_x,cwt_y,scales,s.kappa,s.ts,s.d);

    %% Cui's method
    case 'Cui'
      %% only for morlet wavelet
      %get period from scale
      sinv=1./(scales');
      n1 = length(cwt_x);
      sWxx=smoothwavelet(sinv(:,ones(1,n1)).*(abs(cwt_x).^2),s.ts,periods,s.n_v,scales);
      sWyy=smoothwavelet(sinv(:,ones(1,n1)).*(abs(cwt_y).^2),s.ts,periods,s.n_v,scales);
      Wxy=cwt_x.*conj(cwt_y);
      sWxy=smoothwavelet(sinv(:,ones(1,n1)).*Wxy,s.ts,periods,s.n_v,scales);
   %% method in Matlab wavelet toolbox, reference wcoh.m
    case 'Matlab'
        n1 = length(cwt_x);
        ns = 12;
        invscales = 1./scales;
        invscales = repmat(invscales,1,n1);
        sWxx = smoothCFS(invscales.*abs(cwt_x).^2,scales,s.ts,ns);
        sWyy = smoothCFS(invscales.*abs(cwt_y).^2,scales,s.ts,ns);
        crossCFS = cwt_x.*conj(cwt_y);
        crossCFS = smoothCFS(invscales.*crossCFS,scales,s.ts,ns);
        sWxy = crossCFS;

    %% Aguitar's method
    case 'Aguiar'

     Wxy = cwt_x.*conj(cwt_y); %Formula (23) in [1]

     % Smoothing process
     [n_scales,n_times]= size(Wxy);
     sWxy = zeros(n_scales,n_times);

     %smoothing performed on wavelet power spectrum (23, 26) [1]
     sWxx = abs(cwt_x).^2;
     sWyy = abs(cwt_y).^2;

%%  smoothing in scale
      if any(ismember(supported_windows, Args.window_scale))
          ws_size = fix(s.ws_size/(2*s.n_v));
          if ws_size < 5
              ws_size = 5; % Minimum size rearding Hamming window
          end
          window_function = str2func(Args.window_scale);
          window = window_function(ws_size);
          window = window./(sum(window));
          for i_time = 1:n_times
              sWxx(:,i_time) = conv(sWxx(:,i_time),window,'same');
              sWyy(:,i_time) = conv(sWyy(:,i_time),window,'same');
              sWxy(:,i_time) = conv(Wxy(:,i_time),window,'same');
          end
      else
        warning('WCO: scale smoothing window is ignored');
      end
% % smoothing in time
     if any(ismember(supported_windows, Args.window_time))
         window_function = str2func(Args.window_time);
         for i_scale = 1:n_scales
            % here size considers scale
            ws_size = fix(periods(i_scale)*s.wt_size/s.ts);
            if ws_size < 5
                ws_size = 5; % Minimum size rearding Hamming window
            end
            window = window_function(ws_size);
            window = window./(sum(window));
            sWxx(i_scale,:) = conv(sWxx(i_scale,:),window,'same');
            sWyy(i_scale,:) = conv(sWyy(i_scale,:),window,'same');
            sWxy(i_scale,:) = conv(sWxy(i_scale,:),window,'same');  % reuse values
          end
      else
          warning('WCO: time smoothing window is ignored');
      end

    case 'multiwavelets'
     % check for required multiple cwts
      if order > 1
          sWxy = sum(cwt_x.*conj(cwt_y),3);
          sWxx = sum(abs(cwt_x).^2,3);
          sWyy = sum(abs(cwt_y).^2,3);
      end
    end


 switch Args.coherence
     case 'WCO'
            WCO = (abs(sWxy))./(sqrt(sWxx.*sWyy));
%             WCO(sWxx<1e-10) = 0;  % reset meaningless WCO value to 0
%             WCO(sWyy<1e-10) = 0;  % reset meaningless WCO value to 0
            coherence = WCO;
    case 'cWCO'
            cWCO = sWxy./sqrt(sWxx.* sWyy); % Formula (25) of [1]
%             cWCO(abs(sWxx)<1e-10) = 0;  % reset meaningless WCO value to 0
%             cWCO(abs(sWyy)<1e-10) = 0;  % reset meaningless WCO value to 0
            coherence = cWCO;
    case 'WCO2'
            WCO_2 = (abs(sWxy).^2)./(sWxx.*sWyy);
%             WCO_2(abs(sWxx)<1e-10) = 0;  % reset meaningless WCO value to 0
%             WCO_2(abs(sWyy)<1e-10) = 0;  % reset meaningless WCO value to 0
            coherence = WCO_2;
 end

if (Args.disp_coherence)
    n1 = length(cwt_x);
    time = s.ts:s.ts:s.ts*n1;
    FourierFactor = periods(1)/scales(1);
    figure;
    switch(Args.y_axis)
        case ('scale')
            coi_scale = Args.coi./FourierFactor;
            wave_space_plot(time,Arg.time_series,scales,coherence,...
                           'soi', s.soi,...
                           'coi', coi_scale,...
                           'y_label', Args.y_axis);
        case('period')
            wave_space_plot(time,Args.time_series,periods,coherence,...
                           'soi', s.soi,...
                           'coi', Args.coi,...
                           'y_label', Args.y_axis);
    end
end

% calculate phase locking value and shannon entropy
[phi,E] = PLV(coherence);
varargout={coherence,phi,E};
varargout=varargout(1:nargout);

end
