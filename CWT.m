function varargout = CWT(wx,varargin)
%%  Continuous wavelt transform   
% provides different agirithms to calcualte continious wavelet transformations for real and complex valued signals
  %% REQUIREMENTS
  
    % Continuous Wavelet Transform
    % Generalized Morse Wavelet or Morlet Wavelet
    
  %% INPUTS

    % time_series (length(time_series) by 1 data vector)
    % wavelet type (generalized morse wavelet or morlet wavelet)
    % gamma and beta for generalized morse wavelet
    % frequency limits in herz
    % frequency of interest in herz
    % sampling interval in seconds
    % frequency for mollet wavelet
    % energy density
    % number of octave per wave
    % the order of generalized morse wavelet
    % filter type for gmw
    % normalization type
    % option for display wavelet
    
  %% OUTPUTS
  
    % wavelet spectra
    % scales and period length 
    % cone of inflence 
    % time(in seconds) 
    % morsespace_freq 
    % heisenberg area, radius in time,radius in frequency 
    

  %% REFERENCES
  
    %   [1] Aguiar-Conraria, L. and Soares, M.J. (2010)
    %       "The Continuous Wavelet Transform: A Primer", NIPE Working paper
    %   [2] Lilly, J.M. and Olhede, S. C.(2009), "Higher-Order Properties of
    %       Analytic Wavelets", IEEE Trans. Signal Process. 57 (1), 146-160.
    %   [3]  Lilly, J.M. and Olhede, S. C.(2009), "On the Analytic Wavelet
    %       Transform", eprint arXiv:0711.3834
    %   [4] Torrence, C. and Compo, T.C., "A Prectical Guide to Wavelet
    %	    79(1), 61-78.
    % 	    Analysis" (1998), Bulletin of the American Meteorological Society,

Args=struct('wavelet_type','gmw', ...  % default as gmw, amor for morlet
            'source', 'J',...  % can choose Jlab or Matlab only required for gmw, for morlet it is automatically Matlab
            'filter','nodetrend',....  % option of gmw, can choose detrended or nodetrend to detrended time series before transforming, from Jlab
            'wave_normalization', 'bandpass',...  % options: 'energy', 'bandpass'
                                             ...  % 'energy': unit energy normalization for wavelets uses the unit energy normalization.  The time-
                                             ...  %            domain wavelet energy SUM(ABS(morse wavelet in time domain).^2,1) is then always unity
                                             ...  % 'bandpass': (equal to matlab) the FFT of the wavelet has a peak value of 2 for all frequencies
            'morse_space','ag_manuel',...  % options: 'density', 'length', 'matfreq', 'ag_manuel', 'ag_auto'
                                    ...  % density: uses the the desnity 'energy_density' to calculate the morse space
                                    ...  %         with respect to the 'freq_limits',  from Jlab [REF]
                                    ...  % length: uses length of the time series to calculate the morse space, from Jlab [REF]
                                    ...  % matfreq: mimics the matlab implementation to calculate the morse space
                                    ...  % ag_manuel: inspired by Aguria, uses the manuel specified space
                                    ...  % ag_auto: auto detected space by Agurias
            'y_axis','period',...  % options: 'period', 'scale' used to define the unit of the y achsis in the cwt plot
            'disp_wavelet', false,...  % options: true, false to plot wavelet
            'disp_cwt', false);  % options: true, false to plot the cwt ibcluding coi and soi

 Args=parseargs_special(varargin,Args);

 s = ComplexWaveletTransform; % import the class ComplexWaveletTransform
 n1 = length(wx);
 time = s.ts:s.ts:s.ts*n1;

 switch(Args.wavelet_type)
    case  'gmw'
      % wavelet center frequencies in radians / second to calculate FourierFactor
      [energy_freq ,inst_freq, peak_freq] = gmw_peakf(s);

      % heisenberg - Heisenberg area
      % sd_wavelet_time - radius (standard deviation) in time
      % sd_wavelet_frequency - radius (standard deviation) in the frequency domain
      [heisenberg,sd_wavelet_time,sd_wavelet_frequency,skew] = morsebox(s.ga,s.be);

      % wavelet center frequency in radians / second here set as the peak frequency [1] (49)
      % Matlab uses a different approach to calucale the wavelet center frequency
      %
      wave_cf = peak_freq;  % [rad/s]
      wave_cf_hz = wave_cf/(2*pi*s.ts);  % [Hz]

      FourierFactor = (2*pi)/wave_cf_hz;  % here for [Hz] formular see: [1] (49)
      % Frq = centfrq('morl')./scales;
      switch (Args.source)
          case 'Matlab'
              P_powered = s.be*s.ga; % time-bandwidth product (37) in [2] and Appendix D in [2] or [3]
              % if 1/Args.ts is specified, morsespace_freq of continious wavelet transform is in [Hz]
              % otherwise morsespace_freq is in [cyles / per sample]
              % handling for surrogates
              for i=1:size(wx,2)
                [WX(:,:,i) ,morsespace_freq,coi] = cwt(wx(:,i),1/s.ts, 'morse',...
                                            'VoicesPerOctave',s.n_v,...
                                            'WaveletParameters',[s.ga P_powered],...
                                            'FrequencyLimits',s.freq_limits);
              end

              WXC = WX;
              periods = 1./morsespace_freq;
              scales = periods/FourierFactor;
              coi=1./coi;
              % disp cwt
              if (Args.disp_cwt)
                  figure;
                  cwt(wx(:,1),1/s.ts, 'morse',...
                      'VoicesPerOctave',s.n_v,...
                      'WaveletParameters',[s.ga P_powered],...
                      'FrequencyLimits',s.freq_limits);
              end
        varargout_2={morsespace_freq, heisenberg,sd_wavelet_time,sd_wavelet_frequency};
        otherwise
          switch(Args.morse_space)
          case  'density'
              wrange = s.freq_limits;
              low_freq = wrange(1).*2*pi.*s.ts;  % low freq in [rad / samples] required for morsespace
              high_freq = wrange(2).*2*pi.*s.ts;  % high freq in [rad / samples] required for morsespace
              f_space = morsespace(s.ga,s.be,high_freq,low_freq,s.energy_density)./(s.ts);  % [rad /sample point]
              scales = wave_cf./(f_space.*s.ts);
              morsespace_freq  = wave_cf_hz./scales;  % morsespace in [Hz]
            case 'length'
              f_space = morsespace(s.ga,s.be,n1)./s.ts;
              scales = wave_cf./(f_space.*s.ts);
              morsespace_freq  = wave_cf_hz./scales;  % morsespace in [Hz]
            case 'matfreq'
                frange =  s.freq_limits;  % Obtain the frequency range
                wrange = frange.*s.ts*2*pi;  % Convert frequencies in Hz to radians/sample

                a0 = 2^(1/s.n_v);

                s_min = wave_cf/wrange(2);
                s_max = wave_cf/wrange(1);
                n_o = log2(s_max/s_min);  % number of octaves
                scales = (s_min*a0.^(0:(s.n_v*n_o)))';
                f_space = (wave_cf./(scales.*s.ts))';  % in [in rad / sample points * s]
                morsespace_freq  = wave_cf_hz./scales;  % morsespace in [Hz] == f_space / (2*pi)
            case 'ag_manuel'
                l_samples = n1;
                a0 = 2;  %  [1] (47)
                n_v = s.n_v;  % number of voices per octave: determines the number of samples
                                  % (scales) across this span.
                % own limits
                frange =  s.freq_limits;
                s_min = wave_cf./(frange(2).*(s.ts*2*pi));  % scales here frange is still cyclic (not*2*pi) because FourierFactor is used
                s_max = wave_cf./(frange(1).*(s.ts*2*pi));  % scales here frange is still cyclic (not*2*pi) because FourierFactor is used

                n_o = log2(s_max/s_min);  % number of octaves: determines the span of frequencies (approximate value = 5)

                scale_idx = [0:n_o*n_v];
                scales = (s_min*a0.^(scale_idx/n_v))';
                f_space = wave_cf./(scales.*s.ts);  % in [in rad / sample points * s]
                morsespace_freq  = wave_cf_hz./scales;  % morsespace in [Hz]
           case 'ag_auto'
                l_samples = n1;
                a0 = 2;  % Matlab uses here a value with respect to n_v here [1] (47) is used
                n_v = s.n_v;  % number of voices per octave: determines the number of samples
                                  % (scales) across this span.
                s_min = wave_cf./(2*pi*s.ts);  % smallest scale occurs at the largest frequency of the wavelet here in [rad/s]
                s_max = 2*floor(n_v*log2(l_samples / (2*sd_wavelet_time*s_min)));

                n_o = log2(s_max/s_min);  % number of octaves: determines the span of frequencies (approximate value = 5)

                scale_idx = [0:n_o*n_v];
                scales = (s_min*a0.^(scale_idx/n_v))';
                f_space = wave_cf./(scales.*s.ts);  % in [in rad / sample points * s]
                morsespace_freq  = wave_cf_hz./scales;  % morsespace in [Hz]
          end

          % returns the wavelet in the time domain (w_morse_t)
          % for frequency domain [w_morse_t, w_morse_f ] required
          % f_space are the frequencies at which the Fourier transform of the wavelets reach their maximum amplitudes
           w_morse_t=morsewave(n1,s.K,s.ga,s.be,f_space.*(s.ts),Args.wave_normalization); % for normalization see [2] - [1]
           % to test
           % w_morse_t=morsewave(n1,s.K,s.ga,s.be,wave_cf.*s.ts,s.wave_normalization);

           % wavelet transformation for real compex & real valued signals
           [WX, WXC] = wavetrans(wx, conj(wx), w_morse_t, Args.filter);
           % wavelet transformation for real valued signals
           %WX = wavetrans(wx, w_morse_t, Args.filter);

           %% cone of influence
           % uses radius in time non-dimensionalized with respect to the (radian) peak frequency

           % Fourierfactor hae to be checked
           coiScalar = wave_cf./sd_wavelet_time; % [1] coi scalar
           samples = createCoiIndices(n1);
           coitmp = coiScalar*s.ts*samples;
           coi = 1./coitmp;
           coi(coi>(1/s.ts(1)/2)) = max(f_space);
           coi=1./coi;
           coi = coi';
           %static coi
           %static_coi = coiScalar*Args.ts*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];

           periods = 1./(morsespace_freq);  % in [1/s] not per sample
           %periods = 2*pi./f_space; %  have to be checked! here cycle  2*pi./f_space;
           %scales = periods./FourierFactor;

           % display wavelet
           if (Args.disp_wavelet)
               figure;
               [~, index] = min(abs(f_space-peak_freq));
               plot(real(WX(:,index)));
               hold on
               plot(imag(WX(:,index)));
           end


          varargout_2={morsespace_freq, heisenberg,sd_wavelet_time,sd_wavelet_frequency};
        end
    case  'amor'
           for i=1:size(wx,2)
                 [WX(:,:,i),morsespace_freq,coi]= cwt(wx(:,i),1/s.ts,'amor','VoicesPerOctave',s.n_v,...
                                    'FrequencyLimits',s.freq_limits);
           end

        WXC = WX;
        FourierFactor = (2*pi)/6; % TODO!
        periods = 1./morsespace_freq;
        scales = periods/FourierFactor;
        coi=1./coi;
        varargout_2 = {};
 end
 % ensure [scale time] output
 if size((WX),1) > size((WX),2)
     switch (size(size(WX),2))
         case(4)  %multiwavelets and several time series
          WX = permute(WX, [2 1 3 4]);
           WXC = permute(WXC, [2 1 3 4]);
         case(3)  % multiwavelet or several time series
          WX = permute(WX, [2 1 3]);
           WXC = permute(WXC, [2 1 3]);
         case(2)  % single time series without multiwavelets
             WX = WX';
            WXC = WXC';
         otherwise
             warning('CWT: unknown coherence dimensions');
     end
 end
 if (Args.disp_cwt)
     %power = (abs(WX)).^2 ;
     figure;
     switch(Args.y_axis)
         case ('scale')
             coi_scale = coi./FourierFactor;
             soi = s.soi./FourierFactor;
             wave_space_plot(time,wx,scales,WXC(:,:,1),...
                            'soi', soi,...
                            'coi', coi_scale,...
                            'y_label', Args.y_axis);
         case('period')
             wave_space_plot(time,wx,periods,WXC(:,:,1),...
                            'soi', s.soi,...
                            'coi', coi,...
                            'y_label', Args.y_axis);
     end
 end
 varargout={WX,scales,periods,coi,time, varargout_2};
 varargout=varargout(1:nargout);
end
