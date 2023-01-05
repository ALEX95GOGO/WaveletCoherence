function [Surr,Errors]=IAAWT(x,varargin)

%%  Iterative Amplitude Adjusted Wavelet Transform (IAAWT)

%% REQUIREMENTS

  % dual tree complex wavelet transform
  % Instructions on how to access it are available at
  % http://sigproc.eng.cam.ac.uk/Main/NGK

%% INPUT

  % x (length(x) by 1 data vector)
  % This is a dyadic wavelet transform. If you pass it data that is not of
  % length 2^numlevels where numlevels is an integer given by
  % log(length(x)) / log(2) then it will zero pad up to the next
  % dyadic scale.

  % nSurr (the number of surrogates you wish to generate)

  % accerror (the acceptable error)
  % Both the error from the wavelet part and from the restoration of the
  % original values must be below this value for the code to terminate

  % error_change (acceptable relative error)
  % If the change in accerror from one iteration to the next is less than
  % accerror / error_change then the code terminates

%% OUTPUT

  % Surr (length(x) by numsurrogtes array)
  % The generated surrogate data

% Errors (nSurr by 1 vector)
% The value for toterror at the termination of the algorithm for this
% surrogate

%% REFERENCES

% Keylock, C.J. 2017. Multifractal surrogate-data generation algorithm that preserves pointwise
% Hlder regularity structure, with initial applications to turbulence, Physical Review E 95, 032123,
% https://doi.org/10.1103/PhysRevE.95.032123.


Args=struct('nSurr', 1,...  % the number of surrogates
            'error_change', 100, ...  % acceptable relative error
            'accerror', 0.001);  % acceptable error

Args=parseargs_special(varargin,Args);

%randomisation using Kingsbury wavelets
sizer=size(x);
if sizer(1)==1
    x=x';
    sizer=size(x);
end
exactlevels=log(length(x))/log(2);
numlevels=(floor(log(length(x))/log(2)));

count=1;
if abs(numlevels-exactlevels)>eps
    %We need to zero pad
    temp=zeros(2^(numlevels+1),1);
    count=2^(numlevels+1)-length(x)+1;
    temp(count:length(temp),1)=x;
    x=temp;
end

%wavelet decomp
[Yl,Yh] = dtwavexfm(x,numlevels,'near_sym_b','qshift_b');

% real and amplitudes
for loop1=1:numlevels;
    ampYh{loop1}=abs(Yh{loop1});
    phaseYh{loop1}=angle(Yh{loop1});
end

sortval=sort(x(count:length(x)));
stdval=std(sortval);
num2sort=length(sortval);

for surrloop=1:Args.nSurr
    %make a random dataset and take its imag and phases
    [dummy,shuffind]=sort(rand(num2sort,1));
    shuffind=shuffind-1+count;
    z(shuffind,1) = sortval;
    z(1:count-1)=0;
    [Zl,Zh] = dtwavexfm(z,numlevels,'near_sym_b','qshift_b');

    for loop1=1:numlevels
        newphase{loop1}=angle(Zh{loop1});
    end


    amperror(1)=100;
    waverror(1)=100;
    counter=1;

    while (amperror(counter) > Args.accerror) & (waverror(counter) > Args.accerror)
        %wavelet construction
        oldz=z;
        for loop1=1:numlevels
            newZh{loop1}=ampYh{loop1}.*exp(i.*newphase{loop1});
        end

        z=dtwaveifm(Yl,newZh,'near_sym_b','qshift_b');
        wavdiff=mean(mean(abs(real(z)-real(oldz))));
        waverror(counter+1) = wavdiff/stdval;

        %impose original values
        oldz=z;
        data2sort=z(count:length(x));
        [dummy,shuffind]=sort(real(data2sort));
        shuffind=shuffind-1+count;
        z(shuffind,1) = sortval;
        z(1:count-1)=0;
        ampdiff=mean(mean(abs(real(z)-real(oldz))));
        amperror(counter+1) = ampdiff/stdval;

        %Wavelet step
        [nZl,nZh] = dtwavexfm(z,numlevels,'near_sym_b','qshift_b');
        %get phases and imag
        for loop1=1:numlevels;
            newphase{loop1}=angle(nZh{loop1});
        end

        toterror=amperror(counter+1)+waverror(counter+1);
        oldtoterr=amperror(counter)+waverror(counter);
        if abs((oldtoterr-toterror)/toterror) < (Args.accerror/Args.error_change);
            amperror(counter+1)=-1;
            waverror(counter+1)=-1;
        end
        counter=counter+1;
        clear nZh nZl
    end

    clear amperror waverror
    Surr(:,surrloop)=z(count:length(z));
    Errors(surrloop,1)=toterror;
end
