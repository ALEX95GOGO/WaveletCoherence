function varargout=SigniWCO(inputWCO, cwt_surr_x1, cwt_surr_x2, scales,periods,varargin)

  %% Efficient Calculation of the significant wavelet coherence matrix
  %% REQUIREMENTS

    % WCO

  %% INPUT

    % inputWCO: wavelet coherence matrix of any kind
    % cwt_surr_x1: matrices of wavelet coefficients from cwt of nSur surrogates [w_f_1 w_t_1 iSur_1]
    % cwt_surr_x2: matrices of wavelet coefficients from cwt of nSur surrogates [w_f_2 w_t_2 iSur_2]
    % scales: of inputWCO required for scalogram
    % periods: of inputWCO required for plot
    % see WCO.m for more information

  %% OUTPUT

    % SigniWCO: is the original inputWCO where non significant coherence is set as zero
    % P_WCO is matrix od doubles with equal dimensions to signiWCO, representing the p-values
    % Classifier: matrix of integers, represents the predicted significant
    %             coherence
    % nSurr: is a positive integer representing the number of surrogate wco
    %       calculated from the permuation of each cwt_surr_x1 and cwt_surr_x2

Args=struct('signilevel', 0.01,...  % is the level of significane that determince significance of each
                               ...  % wco coeffient from inputWCO compaired to the corresponding coeffienct of the surrogate wco
            'y_axis', 'period',...  % unit of y axis: scale or period
            'time_series',[],...  % time series recommended for plot
            'coi', [],... % is the coen of influence from inputWCO required for plotWAVE
            'disp_coherence',false,...  % show significant wavelet coherence
            'smoothing','Aguiar',...  % default as cohen, can choose Aguiar CUI or multiwavelets
            'coherence','cWCO',...  % default as cWCO, can choose WCO or WCO2
            'soi',[2.02 12.8],...  % space of interest border (here for period length), unit must be the same as graph
            'kappa', 20,...  % smoothing parameter for cohen
            'd', 1,...  % morlet wavelet shape diameter for cohen
            'dj',1/12,...  % this will do 10 sub-octaves per octave
            'ts',0.1,...  % sampling time
            'disp_wco',false,... % option to plot the coherence,
            'connectivity_vector',[],...  % vector representing the time of expected synchrony
            'ws_size',5,...  % to determine the size of the Hamming window window used for smoothing in the scale direction; size depends on dj and is  given by: ws_size/(2*dj) (with a minimum value of 5)
            'wt_size',5);  % to determine the size of the Hamming window used for smoothing
                            % in the time direction; sizevaries with scale s and
                            % is given by:  wt_size*s/dt (with a minimum value of 5)

Args=parseargs_special(varargin,Args);

%% reduce the required number of surrogates by combining all possible combinations
% of surrogates of cwt_surr_x1 and cwt_surr_x2
cwt_combSurrIdx = combvec([1:size(cwt_surr_x1,3)], [1:size(cwt_surr_x1,3)]); % combntns for equal number of surrogates

%% number of permutated surrogates
nSurr = size(cwt_combSurrIdx,2);

%% compaire coefficients
pWCO=zeros(size(cwt_surr_x1,1), size(cwt_surr_x1,2));

for iSur=1:nSurr
          % ensure [scale time] output
    switch (size(size(cwt_surr_x1),2))
        case(4)  %multiwavelets and several time series
            cwt_x1 = cwt_surr_x1(:,:,cwt_combSurrIdx(1,iSur), :);
            cwt_x2 = cwt_surr_x2(:,:,cwt_combSurrIdx(2,iSur), :);

            % reduce dimensions --> regarding wco for surrogates & multiwavelets
            cwt_x1 = squeeze(cwt_x1);
            cwt_x2 = squeeze(cwt_x2);
        case(3)  % multiwavelet or several time series
            cwt_x1 = cwt_surr_x1(:,:,cwt_combSurrIdx(1,iSur));
            cwt_x2 = cwt_surr_x2(:,:,cwt_combSurrIdx(2,iSur));

            % reduce dimensions --> regarding wco for surrogates & multiwavelets
            cwt_x1 = squeeze(cwt_x1);
            cwt_x2 = squeeze(cwt_x2);
        case(2)  % number of surrogates == 1
            cwt_x1 = cwt_surr_x1(:,:);
            cwt_x2 = cwt_surr_x2(:,:);

        otherwise
                  warning('SigniWCO: unknown coherence dimensions');
    end

    WCOSur = WCO(cwt_x1, cwt_x2, scales,periods,...
                'smoothing',Args.smoothing, ...
                'coherence',Args.coherence, ...
                'soi',Args.soi,...
                'kappa', Args.kappa ,...
                'd', Args.d ,...
                'dj',Args.dj,...
                'ts',Args.ts,...
                'disp_coherence', false, ...
                'ws_size',Args.ws_size, ...
                'wt_size',Args.wt_size);
    pWCO = pWCO + (abs(WCOSur) >= abs(inputWCO));
end

P_WCO = pWCO/nSurr;
Signi_WCO = inputWCO;
Signi_WCO(P_WCO >= Args.signilevel) = 0;
Classifier = P_WCO;
Classifier(P_WCO >= Args.signilevel) = 0;
Classifier(P_WCO < Args.signilevel) = 1;

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
                           'pWCO', P_WCO,...
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
                           'pWCO', P_WCO,...
                           'signilevel', Args.signilevel,...
                           'connectivity_area', connectivity_area,...
                           'y_label', Args.y_axis);
    end
end
varargout={Signi_WCO, P_WCO, Classifier, nSurr};
end
