%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fit an exponential curve to a photoacoustic signal
%   Modified from  "Fit_ua_Spec_measurements.m" by 
%   O. Ogunlade and T.J. Allen and J. Brunker
%   Last modification 19 April 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T = Time scale (e.g. 1e-9 for ns and 1e-6 for microsecond)    
% c = Speed of sound (e.g. 5660 for OD filters, 1485 water )
% Start_point_1 = start point for area of interest for the curve fitting
% End_point_1 = end point for area of interest for the curve fitting
% Start_point_2 = start point for area of interest for the offset
% End_point_2 = end point for area of interest for the offset
% Can show up to 24 fits on the plot

% Plot the data at each wavelength to select start and end points


%clc
clear all %need to clear previous arrays
close all


file = input('Input name of text file (including .txt extension):\n','s');

fid = fopen([file],'r');
fid2 = fopen([file],'r');

%%%%%%%%%%%%%%%%%% Read calibration file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calibrate = menu('Calibrate?','Yes: select calibration file','No');
if calibrate == 1
    [calfile,calpath] = uigetfile('*.txt');
    addpath(calpath);
    calData = dlmread(calfile);
    calFactor = calData(:,2);
end

%%%%%%%%%%%%%%%%%%% Read data file until it reaches time unit  %%%%%%%%%%%%%%%%%

while 1
    tline = fgetl(fid);
    itistimeline = strfind(tline, '#TimeUnit:');
    if isempty(itistimeline) == 0;
        break;
    end
end

TimeUnit = fscanf(fid,'%s',1);

timescale = menu({'Select time scale',...
    ['(file specifies "',TimeUnit,'")']},...
    'nanoseconds (1e-9)','microseconds (1e-6)','Other');
if timescale == 1
    T = 1e-9; TimeUnit = 'ns';
elseif timescale == 2
    T = 1e-6; TimeUnit = '\mus';
elseif timescale == 3
    T = input ('Enter time scale (e.g 1e-9 for ns):\n');
    TimeUnit = input...
        ('Enter string for time scale (e.g. \mus for microsec):\n','s');
end

soundSpeed = menu('Select sound speed','1485 m/s (water)',...
    '5660 m/s (OD filters)','Other');
if soundSpeed == 1
    speed_s = 1485;
elseif soundSpeed == 2
    speed_s = 5660;
elseif soundSpeed == 3
    speed_s = input('Enter sound speed in m/s:\n');
end

datainvert = menu('Invert data?','Yes','No');

%%%%%%%%%%%%%%%%%%% Read file until it reaches data  %%%%%%%%%%%%%%%%%%%%%%

while 1
    tline = fgetl(fid);
    itisdataline = strfind(tline, '#Data:');
    if isempty(itisdataline) == 0;
        break;
    end
end


%%%%%%%%%%%%%%%%%%% Define record length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=0;
while 1
    tline = fgetl(fid2);
    itiswavelengthline = strfind(tline, '#Wavelength:');
    if isempty(itiswavelengthline) == 0;
        break;
    end
    a=a+1;
end
b=a;
while 1
    tline = fgetl(fid2);
    itiswavelengthline = strfind(tline, '#Wavelength:');
    if isempty(itiswavelengthline) == 0;
        break;
    end
    a=a+1;
end

data_pt = 1;
measrows = 2;
meascols = a-b-3;


%%%%%% Read wavelength, pulse energy and data for each wavelength  %%%%%%%%

while ~feof(fid)
    
    mystruct(data_pt).wavelength = fscanf(fid, '%*s \n %g', 1);
    mystruct(data_pt).pdval = fscanf(fid, '%*s %*s %*s \n %g', 1);
    
    % fscanf fills the array in column order, so transpose the results
    mystruct(data_pt).meas  =  fscanf(fid, '%g', [measrows, meascols])';
    data_pt = data_pt + 1;
end
fclose(fid);
fclose(fid2);

%%% Plot the PA signal for each wavelength, option to subtract offset %%%%%

plotrange = menu('Select plot range',...
    'all wavelengths','single wavelength','wavelength range');
if plotrange == 1
    startval = 1; 
    endval = length(mystruct)-1;
    rangetext = 'all-wavelengths';
elseif plotrange == 2
    startval = menu('Select wavelength',...
        mystruct(1:length(mystruct)-1).wavelength);
    endval = startval;
    rangetext = [num2str(mystruct(startval).wavelength),'nm'];
elseif plotrange == 3
    startval = menu('Select start wavelength',...
        mystruct(1:length(mystruct)-1).wavelength);
    endval = menu('Select end wavelength',...
        mystruct(1:length(mystruct)-1).wavelength);
    rangetext = [strcat(num2str(mystruct(startval).wavelength),'-',num2str(mystruct(endval).wavelength)),'nm'];
end
    
cc = jet(length(mystruct));
cc = cc(size(cc,1):-1:1,:);

scrsz = get(0,'ScreenSize');
f1 = figure('Position',[0,0,scrsz(3),scrsz(4)]);
rows = 7;
cols = 6;

% begins from longest wavelength
offset_end = 200;
for ii = startval:endval;
    
    % test to see if software failed to save data
    if isempty(mystruct(ii).meas) == 1
        % do nothing, skip this wavelength
    else
        
        Wavelengths(:,ii) = mystruct(ii).wavelength;
        Pd_energies(:,ii) = mystruct(ii).pdval;
        
        PA_time(:,ii) = mystruct(ii).meas(:,1);
        if datainvert == 1
            PA_amp(:,ii) = -mystruct(ii).meas(:,2);
        elseif datainvert == 2
            PA_amp(:,ii) = mystruct(ii).meas(:,2);
        end
        if calibrate == 1
            PA_amp(:,ii) = (PA_amp(:,ii)./Pd_energies(:,ii))./calFactor(ii);
        end
        
    end
    
    cols2 = cols+2;
    subplot(rows,cols2,rows*cols2);%(rows*cols2-cols2+6):(rows*cols2))
    plot(PA_time(:,:),PA_amp(:,:))
    xlabel(['Time (',TimeUnit,')']);
    set(gcf,'DefaultAxesColorOrder',cc)
    
    subplot(rows,cols,ii)
    plot(PA_time(:,ii),PA_amp(:,ii),'color',cc(ii,:),'linewidth',2)
    xlabel(['Time (',TimeUnit,')']);
    title(['Wavelength: ',num2str(Wavelengths(:,ii)),' nm'])
    hold all
    
    ylims = get(gca,'YLim');
    txt = text(PA_time(offset_end),ylims(2)-(ylims(2)-ylims(1))/10,...
        ['end of offset calc = ',num2str(PA_time(offset_end,ii)),' ',...
        TimeUnit]);
    offset_endplot = plot([PA_time(offset_end) PA_time(offset_end)],...
        ylim,'k');
    offsetval = (mean(PA_amp(1:offset_end,ii)));
    PA_amp_offset(:,ii) = PA_amp(:,ii) - offsetval;
    offsetplot = plot(PA_time(:,ii),PA_amp_offset(:,ii),'color',cc(ii,:));
    
    offset_opt = 2;
    offset_opt = menu('Remove zero offset?',['Subtract offset up to ',...
        num2str(PA_time(offset_end,ii)),' ',TimeUnit],...
        'Choose new end limit for offset subtraction',...
        'Do not subtract offset');
    
    while offset_opt == 2;
        [x,y] = ginput(1); [~,offset_end] = min(abs(PA_time(:,ii)-x));
        set(txt,'String',['end of offset calc = ',...
            num2str(PA_time(offset_end,ii)),' ',TimeUnit],'Position',...
            [PA_time(offset_end),ylims(2)-(ylims(2)-ylims(1))/10]);
        
        delete(offset_endplot)
        offset_endplot = plot([PA_time(offset_end) PA_time(offset_end)],...
            ylim,'k');
        
        delete(offsetplot)
        offsetval = (mean(PA_amp(1:offset_end,ii)));
        PA_amp_offset(:,ii) = PA_amp(:,ii) - offsetval;
        offsetplot = plot(PA_time(:,ii),PA_amp_offset(:,ii),...
            'color',cc(ii,:),'linestyle',':');
        
        offset_opt = menu('Remove zero offset?',...
            ['Subtract offset up to ',num2str(PA_time(offset_end,ii)),...
            ' ',TimeUnit],'Choose new end limit for offset subtraction',...
            'Do not subtract offset');
    end
    
    if offset_opt ~= 3
        PA_amp_offset(:,ii) = PA_amp(:,ii) - offsetval;
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%% Measure amplidue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if calibrate == 1
ifit = startval;
while ifit <= endval
    
    figure(f1)
    subplot(rows,cols,ifit)
    
    [maxAmp,maxAmpI] = max(PA_amp_offset(:,ifit));
    [minAmp,minAmpI] = min(PA_amp_offset(:,ifit));
    
    plot_maxAmp = plot([PA_time(1,ifit),PA_time(end,ifit)],...
            [maxAmp,maxAmp],'k');
    plot_minAmp = plot([PA_time(1,ifit),PA_time(end,ifit)],...
            [minAmp,minAmp],'k');
    
    ampSelect = menu('Correct amplitude?','Yes: keep it','No: select max, then min');
    if ampSelect == 1
    elseif ampSelect == 2
        [~,maxAmp] = ginput(1);
        set(plot_maxAmp,'YData',[maxAmp,maxAmp]);
        
        [~,minAmp] = ginput(1);
        set(plot_minAmp,'YData',[minAmp,minAmp]);
    end
    
    amplitude(ifit) = maxAmp-minAmp;
    
    ifit = ifit+1;
    
end
end

fitopts = menu('Continue to fit exponentials?','Yes','No');
if fitopts == 1

%%%%%%%%%%%%%%%%%%%%%%% Fit exponentials: two methods %%%%%%%%%%%%%%%%%%%%%

start_fit = round(length(PA_time)*1/3);
end_fit = round(length(PA_time)*2/3);

f2 = figure('Position',[0,0,scrsz(3),scrsz(4)]);

ifit = startval;
while ifit <= endval
    
    figure(f1)
    subplot(rows,cols,ifit)
    ylims = get(gca,'YLim');
    
    start_txt = text(0,0,'');
    end_txt = text(0,0,'');
    
    timeaxis = menu(['Zoom in on relevant part of plot ',num2str(ifit),'?'],...
        'Keep current x-axis','Customize x-axis');
    if timeaxis == 2
        [x,y] = ginput(1);
        [~,start_axis] = min(abs(PA_time(:,ifit)-x));
        set(start_txt,'Position',[PA_time(start_axis,ifit),...
            ylims(1)+(ylims(2)-ylims(1))/5],'String',...
            ['axis startpoint = ',num2str(PA_time(start_axis,ifit)),' ',...
            TimeUnit]);
        
        [x,y] = ginput(1);
        [~,end_axis] = min(abs(PA_time(:,ifit)-x(1)));
        set(end_txt,'Position',[PA_time(end_axis,ifit),...
            ylims(1)+(ylims(2)-ylims(1))/10],'String',...
            ['axis endpoint = ',num2str(PA_time(end_axis,ifit)),' ',...
            TimeUnit]);
        
        xlim([PA_time(start_axis,ifit) PA_time(end_axis,ifit)]);
        ylims = get(gca,'YLim');
        
        start_fit = start_axis;
        end_fit = end_axis;
    end
    
    
    fit_startplot = plot([PA_time(start_fit,ifit)...
        PA_time(start_fit,ifit)],ylim,'k');
    set(start_txt,'Position',[PA_time(start_fit),...
        ylims(1)+(ylims(2)-ylims(1))/5],'String',...
        ['fit startpoint = ',num2str(PA_time(start_fit,ifit)),' ',...
        TimeUnit]);
    
    fit_endplot = plot([PA_time(end_fit,ifit)...
        PA_time(end_fit,ifit)],ylim,'k');
    set(end_txt,'Position',[PA_time(end_fit),...
        ylims(1)+(ylims(2)-ylims(1))/10],'String',...
        ['fit endpoint = ',num2str(PA_time(end_fit,ifit)),' ',...
        TimeUnit]);
    
    fit_opt = 2;
    fit_opt = menu('Choose fitting region',['Fit data from ',...
        num2str(PA_time(start_fit,ifit)),' to ',...
        num2str(PA_time(end_fit,ifit)),' ',TimeUnit],...
        'Choose new start and end points');
    
    while fit_opt == 2;
        [x,y] = ginput(1);
        [~,start_fit] = min(abs(PA_time(:,ifit)-x));
        set(start_txt,'String',['fit startpoint = ',...
            num2str(PA_time(start_fit,ifit)),' ',TimeUnit],'Position',...
            [PA_time(start_fit,ifit),ylims(1)+(ylims(2)-ylims(1))/5]);
        
        [x,y] = ginput(1);
        [~,end_fit] = min(abs(PA_time(:,ifit)-x(1)));
        set(end_txt,'String',['fit endpoint = ',...
            num2str(PA_time(end_fit,ifit)),' ',TimeUnit],'Position',...
            [PA_time(end_fit,ifit),ylims(1)+(ylims(2)-ylims(1))/10]);
        
        delete(fit_endplot)
        fit_endplot = plot([PA_time(end_fit,ifit)...
            PA_time(end_fit,ifit)],ylim,'k');
        
        delete(fit_startplot)
        fit_startplot = plot([PA_time(start_fit,ifit)...
            PA_time(start_fit,ifit)],ylim,'k');
        
        fit_opt = menu('Choose fitting region',['Fit data from ',...
            num2str(PA_time(start_fit,ifit)),' to ',...
            num2str(PA_time(end_fit,ifit)),' ',TimeUnit],...
            'Choose new start and end points');
    end
    
    fit_time = PA_time(start_fit:end_fit,ifit)*T;
    fit_amp =  PA_amp(start_fit:end_fit,ifit);
    
    x = fit_time;
    y = fit_amp;
    
    % -------------------- Fitting method 1 ------------------------------
    % Here we define the exponential
    % Equation of the form y = a*exp(b*x)
    
    [fresult_1, gof] = fit(x,y,'exp1', 'Start', [0,0]);
    conf_1 = confint(fresult_1);
    
    fitted_a(ifit) = fresult_1.a;
    fitted_b(ifit) = fresult_1.b;
    
    fitted_mua_1(ifit) = fitted_b(ifit)./(speed_s*1000); %ua in mm^(-1)
    Lconf_mua_1(ifit) = conf_1(1,2)./(speed_s*1000); %lower 95% confidence limit
    Uconf_mua_1(ifit) = conf_1(2,2)./(speed_s*1000); %upper 95% confidence limit
    
    fitted_a_1_unscaled(ifit) = (fresult_1.a - mean(fit_amp))/std(fit_amp);
    
    
    % store performance metric of the fit method
    fit_quality_1(ifit) = gof.rsquare;
    
    % ---------------------- Fitting method 2 ----------------------------
    % Here we define the exponential
    % Equation of the form a*exp(b*x)+c
    
    
    
    mean_x = mean(x);  std_x = std(x);
    x = (x-mean_x)./std_x;
    
    c_est = 0;
    b_est = (log((y(end)-c_est)/(y(1)-c_est))) / (x(end)-x(1));  % b
    a_est = mean((y-c_est)./exp(b_est*x));
    
    
    
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Start', [a_est b_est c_est],'Algorithm',...
        'Trust-Region','Normalize','on'); % [exp coeff, exp term, offset]
    Myfittype = fittype('FK* exp(Fmua_times_speed_s* (x)) + Offset',...
        'indep','x',...
        'depen','y',...
        'coef',{'FK','Fmua_times_speed_s', 'Offset'},...
        'options',s);
    
    % Produce Fit in "ActualFit"
    [AS,AFK] = fit(fit_time,fit_amp,Myfittype);
    conf_2 = confint(AS); %confidence intervals
    
    % Create data with calculated parameters from ActualFit
    fresult_2 = AS.FK * exp (AS.Fmua_times_speed_s*(x)) + AS.Offset;
    
    
    fitted_b_2(ifit) = AS.FK;
    fitted_c(ifit) = AS.Offset;
    
    
    % since the data has been centered and scaled,
    % the actual sound speed is now
    
    speed_s_scaled = speed_s*std(fit_time) + mean(fit_time);
    fitted_mua_2(ifit) = AS.Fmua_times_speed_s/(speed_s_scaled*1000);
    Lconf_mua_2(ifit) = conf_2(1,2)/(speed_s_scaled*1000); %lower 95% confidence limit
    Uconf_mua_2(ifit) = conf_2(2,2)/(speed_s_scaled*1000); %upper 95% confidence limit
    %ua in mm^(-1)
    
    fitted_a_2_unscaled(ifit) = (AS.FK - mean(fit_amp))/std(fit_amp);
    fitted_c_unscaled(ifit) = (AS.Offset - mean(fit_amp))/std(fit_amp);
    
    
    % store performance metric of the fit method
    fit_quality_2(ifit) = AFK.rsquare;
    
    % plot both original data and both fitted data to see match
    figure(f2)
    subplot(rows,cols,ifit)
    plot(fit_time,fit_amp, 'r-', fit_time,fresult_1(fit_time), ...
        'g-', fit_time,fresult_2, 'b-');
    xlabel(['Time (',TimeUnit,')']);
    if ifit == 0
        leg = legend('Measured data', 'fitted to: a*exp(mua*t)',...
            'fitted to: a*exp(mua*t) + c', 'Location','NorthOutside');
        set(leg,'FontSize',8)
    end
    
    % choose the best of the fit methods for the fitted mua
    if fit_quality_2(ifit)>fit_quality_1(ifit)
        fitted_mua_best(ifit) = fitted_mua_2(ifit);
        Lconf_mua_best(ifit) = Lconf_mua_2(ifit);
        Uconf_mua_best(ifit) = Uconf_mua_2(ifit);
        title([num2str(Wavelengths(:,ifit)),'; a*exp(mua*t) + c'])
    else
        fitted_mua_best(ifit) = fitted_mua_1(ifit);
        Lconf_mua_best(ifit) = Lconf_mua_1(ifit);
        Uconf_mua_best(ifit) = Uconf_mua_1(ifit);
        title([num2str(Wavelengths(:,ifit)),'; a*exp(mua*t)'])
    end
    continueopts = menu(['Repeat plot ',num2str(ifit),'?'],...
        ['Yes, repeat plot ',num2str(ifit)],['No, continue to plot ',num2str(ifit+1)],'End script');
    if continueopts == 1
    elseif continueopts == 2
        ifit = ifit+1;
    elseif continueopts == 3
        endval = ifit;
        break
    end
end

figure(f2)
subplot(rows,cols2,rows*cols2);%(rows*cols2-cols2+6):(rows*cols2))
errorbar(Wavelengths(startval:endval),fitted_mua_best(startval:endval),...
    fitted_mua_best(startval:endval)-Lconf_mua_best(startval:endval),...
    Uconf_mua_best(startval:endval)-fitted_mua_best(startval:endval),'Marker','o')
legend(['Best mua fit (sound speed: ',num2str(speed_s),' m/s)'],'Location','best')

saveas(f2,strcat(file(1:length(file)-4),'_',rangetext,'.fig'))

fitted_mua_percm = fitted_mua_best(startval:endval);
error_bar_length = Uconf_mua_best(startval:endval)-Lconf_mua_best(startval:endval);

end

wavelengths_nm = Wavelengths(startval:endval);
amplitudes = amplitude(startval:endval);

if calibrate == 1
    subplot(rows,cols2,rows*cols2);%(rows*cols2-cols2+6):(rows*cols2))
    hold on
    plot(Wavelengths(startval:endval),amplitude(startval:endval),'r','Marker','x')
    
    legend(['Best mua fit (sound speed: ',num2str(speed_s),' m/s)'],'Amplitude','Location','best')
    
    saveas(gcf,strcat(file(1:length(file)-4),'_',rangetext,'_amplitudes.fig'))
    
    if fitopts == 1
        T = table(wavelengths_nm,fitted_mua_percm,error_bar_length,amplitudes);
        writetable(T,strcat(file(1:length(file)-4),'_',rangetext,'_fits,amplitudes.txt'));
        writetable(T,strcat(file(1:length(file)-4),'_',rangetext,'_fits,amplitudes.csv'));
    else
        wavelengths_nm = wavelengths_nm';
        amplitudes = amplitudes';
        T = table(wavelengths_nm,amplitudes);
        writetable(T,strcat(file(1:length(file)-4),'_',rangetext,'_amplitudes.txt'));
        writetable(T,strcat(file(1:length(file)-4),'_',rangetext,'_amplitudes.csv'));
    end
else
    T = table(wavelengths_nm,fitted_mua_percm,error_bar_length);
    writetable(T,strcat(file(1:length(file)-4),'_',rangetext,'_fits.txt'));
    writetable(T,strcat(file(1:length(file)-4),'_',rangetext,'_fits.csv'));
end




