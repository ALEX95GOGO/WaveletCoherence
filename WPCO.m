function varargout = WPCO(cwt_1,cwt_2)
%% Wavelet Phase Coherence

%% INPUTS

  % cwt_1 wavelet transform of signal 1
  % cwt_2 wavelet transform of signal 2

%% OUTPUTS

  % wavelet phase coherence

%% REFERENCE
    % Bloomfield et al. Wavelet_Phase_Coherence_Analysis__Application_to_a_Quiet_Sun_Magnetic_Element
    % Sheppard et al. Testing for time-localized coherence in bivariate data

% required for different sizes
f_min = min([size(cwt_1,1),size(cwt_2,1)]);
cwt_1=cwt_1(1:f_min,:);
cwt_2=cwt_2(1:f_min,:);

% angle
phi1=angle(cwt_1);
phi2=angle(cwt_2);
exp_phase_difference = exp(1i*(phi1-phi2));

WPCO=zeros(1,f_min)*NaN;
phase_difference=zeros(1,f_min)*NaN;
for f=1:f_min
    cphexp=exp_phase_difference(f,:);
    cphexp=cphexp(~isnan(cphexp));
    no_lack=length(find(cwt_1(f,:)==0 & cwt_2(f,:)==0));
    contains_lack=length(cphexp);
    if contains_lack>0
        phph=mean(cphexp)-no_lack/contains_lack;
        WPCO(f)=abs(phph);
        phase_difference(f)=angle(phph);
    end
end

varargout={WPCO',phase_difference'};
varargout=varargout(1:nargout);

end
