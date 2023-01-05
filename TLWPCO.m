function TPC = TLWPCO(cwt_1,cwt_2,scales, periods, varargin)
%% Time-localized Wavelet Phase Coherence

%% INPUTS

  % cwt_1 wavelet transform of signal 1
  % cwt_2 wavelet transform of signal 2

%% OUTPUTS

  % TPC wavelet phase coherence

%% REFERENCE

  % Sheppard et al. Testing for time-localized coherence in bivariate data

  
 Args = struct( 'coi',[],...  % coen of influence recommended for plot
                'tx',[],...  % time series recommended for plot
                'ty',[],...  % time series recommended for plot'ts', 0.1,...  % sampling rate in [s]
                'disp_coherence', false, ... %option to plot the coherence
                'y_axis','period',... % can choose Scale and Period
                'plotsmooth', 0.5, ...  % smoothing factor of plot [0 : 1]
                'FOI_h', 2.02, ...  % frequency of interest high border (here for period length), unit must be the same as graph
                'FOI_l', 12.8, ...  % frequency of interest lower border
                'ts',0.1,...  % sampling time [s]
                'n_cycles', 10);  % number of cycles for calculating TPC (determines adaptive
                                % window length
Args=parseargs_special(varargin,Args);

[n_f,n_t]=size(cwt_1);

cross_phase_coherence = exp(1i*angle(cwt_1.*conj(cwt_2)));
zero_cph = cross_phase_coherence;
zero_cph(isnan(zero_cph))=0; 
cumm_cph = [zeros(n_f,1),cumsum(zero_cph,2)];
TPC = zeros(n_f,n_t)*NaN;

for fn=1:n_f
    cross_phase_current = cross_phase_coherence(fn,:); 
    cumm_cph_current = cumm_cph(fn,:);
    tn1 = find(~isnan(cross_phase_current),1,'first'); 
    tn2 = find(~isnan(cross_phase_current),1,'last');

    window = round((Args.n_cycles*periods(fn))./Args.ts);  % attention here periods are assumed to be in [s] not radiants
    window = window+1-mod(window,2);
    hw =  floor(window/2);

    if ~isempty(tn1+tn2) && window<=tn2-tn1
        local_phase_coherence=abs(cumm_cph_current(tn1+window:tn2+1)-...
                              cumm_cph_current(tn1:tn2-window+1))/window;
        TPC(fn,tn1+hw:tn2-hw)=local_phase_coherence;
    end
end
if (Args.disp_coherence)
   time = Args.ts:Args.ts:Args.ts*n_t;
  if (size(Args.tx) & size(time') & size(Args.ty) == size(time') & size(Args.coi) == size(time'))
    figure;
       
    % consider unit
    FourierFactor = periods(1)/scales(1);
    FOI_h = time;
    FOI_l = time;

    switch (Args.y_axis)
           case('scale')
               wavespecplot_s(time,Args.tx,Args.ty,scales,TPC', Args.plotsmooth);
               ylabel('scale');
               xlabel('time (s)');
               hold on
               % cone-of-influence, anything "below" is dubious
               %plot(time,Args.coi./FourierFactor,'w',...
               %       'LineWidth',2);
                baselevel = max(Args.coi);
                coi_area = area(time,Args.coi./FourierFactor,baselevel);
                
                FOI_h(:) = Args.FOI_h./FourierFactor;
                FOI_l(:) = Args.FOI_l./FourierFactor;
                
           case ('period')
                wavespecplot_s(time,Args.tx,Args.ty,periods,TPC', Args.plotsmooth);
                ylabel('period length(s)');
                xlabel('time (s)');
                hold on

                % cone-of-influence, anything "below" is dubious
                %plot(time,Args.coi,'w',...
                %              'LineWidth',2);
                baselevel = max(Args.coi);
                coi_area = area(time,Args.coi,baselevel);
                
                FOI_h(:) = Args.FOI_h;
                FOI_l(:) = Args.FOI_l;
                
    end
 
    coi_area.FaceAlpha = 0.6;
    coi_area.FaceColor = 'w';
    coi_area.LineStyle = 'none';
  
    plot(time,FOI_h,'w',...
                  'LineWidth',1);
    plot(time,FOI_l,'w',...
                  'LineWidth',1);
  else
    warning('PLOT Coherence, time series or coi is missing or missmatched');
  end
end
