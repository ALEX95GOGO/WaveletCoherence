function varargout=SigniWCOCutOff(inputWCO, cut_off, scales, periods,varargin)

  %% Efficient Calculation of the significant wavelet coherence matrix with a cut off value
  %% REQUIREMENTS

    % WCO

  %% INPUT

    % inputWCO: wavelet coherence matrix of any kind
    % cut_off: double, determines which coherence is significant
    % scales: of inputWCO required for scalogram
    % periods: of inputWCO required for plot
    % see WCO.m for more information

  %% OUTPUT

    % SigniWCO: is the original inputWCO where non significant coherence is set as zero
    % P_WCO is matrix od doubles with equal dimensions to signiWCO, representing the p-values
    % Classifier: matrix of integers, represents the predicted significant
    %             coherence

Args=struct('y_axis', 'period',...  % unit of y axis: scale or period
            'time_series',[],...  % time series recommended for plot
            'coi', [],... % is the coen of influence from inputWCO required for plotWAVE
            'disp_coherence',false,...  % show significant wavelet coherence
            'soi',[2.02 12.8],...  % space of interest border (here for period length), unit must be the same as graph
            'ts',0.1,...  % sampling time
            'disp_wco',false,... % option to plot the coherence,
            'connectivity_vector',[]...  % vector representing the time of expected synchrony
            );
Args=parseargs_special(varargin,Args);


Signi_WCO = inputWCO;
Signi_WCO(inputWCO < cut_off) = 0;
Classifier = Signi_WCO;
Classifier(Signi_WCO > 0) = 1;

if (Args.disp_coherence)
    n1 = length(inputWCO);
    time = Args.ts:Args.ts:Args.ts*n1;
    FourierFactor = periods(1)/scales(1);
    figure;
    switch(Args.y_axis)
        case ('scale')
            coi_scale = Args.coi./FourierFactor;
            %% create connectivity area from vector and set only soi as relevant
            connectivity_area = zeros(size(inputWCO,1), size(inputWCO,2));
            idx = find(scale >= min(Args.soi) & scale <= max(Args.soi));
            connectivity_area(idx,:) = repmat(Args.connectivity_vector', length(idx),1);
            %% plot
            wave_space_plot(time,Args.time_series,scales,inputWCO,...
                           'soi', Args.soi,...
                           'coi', coi_scale,...
                           'pWCO', Classifier,...
                           'signilevel', Args.signilevel,...
                           'connectivity_area', Args.connectivity_area,...
                           'y_label', Args.y_axis);
         case('period')
            %% create connectivity area from vector and set only soi as relevant
            connectivity_area = zeros(size(inputWCO,1), size(inputWCO,2));
            idx = find(periods >= min(Args.soi) & periods <= max(Args.soi));
            connectivity_area(idx,:) = repmat(Args.connectivity_vector', length(idx),1);
            %% plot
            wave_space_plot(time,Args.time_series,periods,inputWCO,...
                           'soi', Args.soi,...
                           'coi', Args.coi,...
                           'pWCO', Classifier,...
                           'connectivity_area', connectivity_area,...
                           'y_label', Args.y_axis);
    end
end
varargout={Signi_WCO, Classifier};
end
