function wave_space_plot(time, time_series, frequency_space, WXC, varargin)
%WAVESPECPLOT  Plot of wavelet spectra together with time series.

%% INPUT
    % time
    % time_series a vector or matrix of time series
    % frequency_space vector of doubles can be in Hz [frequency], periods [s] or scale [s]
    % WXC wavelet transform regarding frequency_space and time_series


  Args=struct('coi',[],...  % coen of influence
              'soi',[],...  % space of interest
              'pWCO',[],...  % representing the p-values
              'signilevel', 0.01,...  % Levels of significance 5%
              'connectivity_area', [],...  % Vectors representing the area of expected synchrony
              'y_label','scale [rad]');  % label of wave space subplot
  Args=parseargs_special(varargin,Args);

  time_min=min(time);
  time_max=max(time);

  if isreal(time_series)
    lim_y = maxmax(abs(time_series))*1.1;
  else
    lim_y = max([maxmax(abs(real(time_series))) maxmax(abs(imag(time_series)))])*1.1;
  end


  subplot(2,1,1)
%   if ~isreal(time_series)
%     sub_p1 = plot(ii,real(time_series)); hold on
%     sub_p2 = plot(ii,imag(time_series));
%     sub_plot = [sub_p1;sub_p2];
%   else
    sub_plot = plot(time,time_series);
 % end

  xlim([time_min,time_max])
  axis 'auto x'
  ylim([-lim_y lim_y])
  axis 'auto y'
  ylabel(gca,strcat(char(8710), 'Concentration [(mmMol/l)*mm]'),...
        'FontName','calibri','FontSize',9);
  grid off

  subplot(2,1,2)
  if ~isreal(WXC)
      WXC=abs(WXC);
  end

  pcolor(time,frequency_space,WXC),shading interp

  xlim([time_min,time_max])
  xlabel(gca,'Time (seconds) ','FontName','calibri','FontSize',9);
  ylabel(gca,Args.y_label,...
        'FontName','calibri','FontSize',9);
  hold on
  set(gca,'YDir','reverse');
  set(gca,'tickdir','out')
  set(gca,'yscale','log');
  colorbar;

%% cone of influence
  if length(Args.coi) == length(time)
    baselevel = max(Args.coi);
    coi_area = area(time,Args.coi,baselevel);
    coi_area.FaceAlpha = 0.6;
    coi_area.FaceColor = 'w';
    coi_area.LineStyle = 'none';
  end
%% scale of interest
  if length(Args.soi) > 0
    soi_l = time;
    soi_h = time;
    soi_l(:) = min(Args.soi);
    soi_h(:) = max(Args.soi);
    plot(time,soi_l,'w',...
                  'LineWidth',1);
    plot(time,soi_h,'w',...
                  'LineWidth',1);
  end
  %% connectivity area
  if length(Args.connectivity_area) > 0
%       colormap = gray(2);
%     plot_connectivity_area = pcolor(time,frequency_space,Args.connectivity_area);
%     plot_connectivity_area.EdgeAlpha = 0;
%     plot_connectivity_area.FaceAlpha = 0.15;
    contour(time,frequency_space,Args.connectivity_area,[1, 1], 'magenta','LineWidth',1);
  end
  %% significance
  if sum(size(Args.pWCO)) > 0
      contour(time,frequency_space,Args.pWCO,[Args.signilevel, Args.signilevel],...
     'k-','LineWidth',0.8);
  end
hold off;
