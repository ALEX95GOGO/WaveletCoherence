%% Local Injected Synchrony Simulation (LISS)

%% INPUT

  % ts_x:, ts_y vector of real valued time series representing the signals
  % ts: sampling time [s]

%% OUTPUT

  % ts_X: vector of doubles representing the original time series ts_x
  % ts_Y: vector of doubles representing the time series with injected loclized signal components of ts_x
  % Labes: vector of integers representing the defined "ground truth"
  % Labels_Shifted: vector of integers representing the shifted window
  % RECT: vector integers the locations of artifical synchrony

%% REQUIREMENTS

  % WhiteGaussianNoise

%% REFERENCE

function varargout = LISS(ts_x, ts_y, varargin)

Args=struct('disp_info', false,...  % to enable periodogram and wco           
            'filter', 'bandstop');  % options: knockout, bandstop, otherwise
                                    % knowckout uses bandpass reduce irrelevant FOI parts of manipulated signal
                                    % bandstop with zero-phase shift will force FOI of manipulated signal to zero
                                    % otherwise the FOI is preserved

Args=parseargs_special(varargin,Args);
s = LocalInjectedSynchronySimulation; % import the class LocalInjectedSynchronySimulation

%% simple length adjustment
if length(ts_x) ~= length(ts_y)
  lmin = min(length(ts_x), length(ts_y));
  ts_x = ts_x(1:lmin);
  ts_y = ts_y(1:lmin);
end
l = length(ts_x);

%% use RECT for partial identical signals
s.tr = s.tr';
s.T = s.T';
RECT = zeros(l, 1);
RECT_seperated =  zeros(l, length(s.T));
for i = 1:length(s.T)
    RECT(s.tr(i) - s.T(i)/2: s.tr(i) + s.T(i)/2) = 1;
    RECT_seperated(s.tr(i) - s.T(i)/2: s.tr(i) + s.T(i)/2, i) = 1;  % creates rect for each window
end

%% get labeled area
% only joint areas of a specific window are labelled as true
% joint areas between of different windows are not true
Label = RECT_seperated;
Label_Shifted = circshift(RECT_seperated, s.shift_time);
Label = Label.*Label_Shifted;  % only joint set is assumed as true
Label = sum(Label,2);

%% filters only specific frequency range, regarding the frequency shift
ts_x_w = bandpass(ts_x,s.FOI,1/s.ts);
ts_y_w = bandpass(ts_y,s.FOI,1/s.ts);

switch Args.filter
  case 'bandstop'
 % while a notch filter notches out only a specific frequency, here we have
 % to filter a frequency
    ts_X = ts_x;

%% self designed filter
%   bandstop = designfilt('bandstopiir','FilterOrder',56, ...  % influences width can be optimized
%          'HalfPowerFrequency1',Args.FOI(1),'HalfPowerFrequency2',Args.FOI(2), ...
%          'DesignMethod','butter',...
%          'SampleRate',1/Args.ts);
%   fvtool(bandstop)
%   ts_y_n = filtfilt(notch_filter, ts_y);
    ts_y_bandstop = bandstop(ts_y,s.FOI,1/s.ts,'ImpulseResponse', 'auto');  % iir or fir filter
    ts_y_foi = RECT.*ts_x_w;
    ts_y_foi = circshift(ts_y_foi, s.shift_time);  % shift signal in time
    ts_Y = ts_y_bandstop + ts_y_foi;  % reconstruct signal
    if Args.disp_info
        figure;
        bandstop(ts_y,s.FOI,1/s.ts,'ImpulseResponse', 'auto')
        title('ts_y_bandstop: manipulated');
     end
   case 'knockout'  % == knockout
    ts_X = ts_x;
    ts_y_foi = RECT.*ts_x_w;
    ts_y_foi = circshift(ts_y_foi, s.shift_time);  % shift signal parts in time

    ts_Y = ts_y-ts_y_w+ts_y_foi;  % reconstruct signal
  otherwise
   ts_X = ts_x;
   ts_y_foi = RECT.*ts_x_w;
   ts_y_foi = circshift(ts_y_foi, s.shift_time);  % shift signal parts in time
   ts_y_foi = ts_y_foi+(1-RECT).*ts_y_w;
   ts_Y = ts_y-ts_y_w+ts_y_foi;  % reconstruct signal, FOI is partly preserved
end

%% add white Gaussian noise to signal
if s.snr > 1
  ts_Y = WhiteGaussianNoise(ts_Y,s.snr);  % adds additive gaussian (white) noise to the manipulated signal with 10dB SNR
end

if Args.disp_info
   figure;
   plot_area = area(RECT);
   plot_area.FaceAlpha = 0.05;
   plot_area.FaceColor ='black';
   hold on;
   plot_area_shifted = area(circshift(RECT, s.shift_time));
   plot_area_shifted.FaceAlpha = 0.05;
   plot_area_shifted.FaceColor ='blue';
   hold on;
   plot (circshift(RECT,s.shift_time));
   hold on;
   plot(ts_X');
   hold on;
   plot(ts_Y');
   ylim([minmin([ts_X, ts_Y])*(1-0.1), maxmax([ts_X, ts_Y])*(1+0.1)]);
   title('marker');
   f = figure;
   p = uipanel('Parent',f,'BorderType','none');
   p.Title = 'Periodogram';
   p.TitlePosition = 'centertop';
   p.FontSize = 12;
   p.FontWeight = 'bold';
   subplot(2,1,1,'Parent',p)
   periodogram(ts_X,hamming(length(ts_y_foi)),256,10)
   title('Periodogram: ts_X');
   subplot(2,1,2,'Parent',p)
   periodogram(ts_Y,hamming(length(ts_y_foi)),256,10)
   title('Periodogram: ts_Y');
end
varargout={ts_X, ts_Y, Label, RECT, Label_Shifted};
varargout=varargout(1:nargout);
end
