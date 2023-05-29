% Alyssa Holman
% Last edited: 10/21/2020
% Calcium Analysis

%% Set paths and open files
% First, set the path below to open the images of interest
clc
clear

%%% CHANGE THIS %%%
% cd('/Users/Alyssa/Desktop/Ca2+/09222020 - mkate gcamp/1900_scan');
cd('');

%%% MAKE A FOLDER ONLY WITH LINESCAN IMAGES %%%
% import your images
imagefiles = dir('*.tif')  % information about each .tif  
nfiles = length(imagefiles);    % number of files found

% get info from files
for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   cropcurrent = currentimage;
   images{ii} = cropcurrent; % extract pixel intensity for each image
end
%% Plot ca2+ flux for each .tif file
% Here, we will take the mean of each row and plot
% and quickly visualize the flux

%%% SET TIME FOR LINESCANNING IN MS %%%
millisec = 30000;
frames = length(mean(images{1},2));
conv = millisec/frames; % gets conversion based on time / frames

for x=1:length(images)
    sub = length(images)/4; % to make a grid of plots in 4 columns
    subplot(ceil(sub),4,x); % to specify number of rows
    y = (1:size(images{x},1)); % x-axis is time in px
    z = mean(images{x},2); % y-axis is intensity in px
    plot(conv*y,z);
    title(imagefiles(x).name);
    xlabel('Time in MS'); % currently, unit is in px not ms
    ylabel('Intensity'); % currently, unit is px intensity not F/Fo
end

%% Smoothen out curve to get better metrics

for x=1:length(images)
    sub = length(images)/4; % to make a grid of plots in 4 columns
    subplot(ceil(sub),4,x); % to specify number of rows
    y = (1:size(images{x},1)); % x-axis is time in px
    m_z = mean(images{x},2); % y-axis is intensity in px
    
    % smoothing function
    % increasing 0.05 makes more smooth
    z = smooth(m_z, 0.05, 'loess'); 
    
    plot(y,z);
    title(imagefiles(x).name);
    xlabel('Time in Pixels'); % currently, unit is in px not ms
    ylabel('Intensity'); % currently, unit is px intensity not F/Fo
end


%% Finding peaks and minima for intensity values
% Here, we want to find the maxima and minima so we can calculate the
% features within the ca2+ flux plots

for x=1:length(images)
    y = 1:size(images{x},1);
    m_z = mean(images{x},2);
    z = smooth(m_z, 0.05, 'loess');
    
    inv = -z;
    
    % preliminary peaks to get ranges
    %%% tune the 'MinPeakDistance' parameter %%%
    [IntMaxs_tmp, TimeMaxs_tmp] = findpeaks(z, y,'MinPeakDistance', 0.5, 'MinPeakProminence', 2);
    [IntMins_temp, TimeMins_tmp] = findpeaks(-z, y,'MinPeakDistance', 0.5, 'MinPeakProminence', 2);
    IntMins_tmp = -IntMins_temp;
    
    % get actual max and min values based on baseline
    %%% INCREASE multiplier if too many peaks, DECREASE if too few %%%
    %%% tune the 'MinPeakDistance' parameter %%%
    prom_min = 0.1*(mean(IntMaxs_tmp) - mean(transpose(IntMins_tmp)));
    prom_max = 0.1*(mean(IntMaxs_tmp) - mean(transpose(IntMins_tmp)));
    
    [IntMaxs, TimeMaxs, width, prom] = findpeaks(z, y,'MinPeakDistance', 0.5, 'MinPeakProminence', prom_max);
    [IntMins_neg, TimeMins] = findpeaks(-z, y,'MinPeakDistance', 0.5, 'MinPeakProminence', prom_min);
    IntMins = -IntMins_neg;
    
    % save values
    IntMax{x} = IntMaxs; % save Intensity Maximum vals (i.e. intensity at peak)
    TimeMax{x} = TimeMaxs; % save Time Maximum vals (i.e. time of peak)
    IntMin{x} = IntMins; % save Intensity Minimum vals (i.e. intensity at minima)
    TimeMin{x} = TimeMins; % save Time Minimum vals (i.e. time of minima)
    Proms{x} = prom;
    width_fin = width*conv;
    Widths{x} = width_fin;
    
    % plot to make sure peaks and minima are correct
    sub = length(images)/4;
    subplot(ceil(sub),4,x); 
    hold on
        plot(y, z)
        title('Resting and Maximum Points');
        plot(TimeMins, IntMins, '*') % adds * where mins are
        plot(TimeMaxs, IntMaxs, '*') % adds * where maxs are
        title(imagefiles(x).name);
        xlabel('Time in Pixels');
        ylabel('Intensity');
    hold off
    
end

%% Finding F and Fo values from intensity values
% First, calculate resting flux and use this to get F/Fo

for x=1:length(images)
    
    % get resting values
    restings = mean(IntMin{x}(2:end,1)); % skipped first val because sometimes not accurate
    resting{x} = restings; % save resting values
    
    % get deltaF values
    m_z = mean(images{x},2);
    z = smooth(m_z, 0.05, 'loess');
    deltF = [];
    
    for i = 1:size(z,1)
        v = (((z(i,1) - restings)/restings)+0.03); % take (intensity - resting) / resting
        deltF = [deltF, v]; % these are the delta(F/Fo) values
    end
    
    % save deltaF values
    deltaF{x} = deltF;
    
    % plot with deltaF values
    sub = length(images)/2;
    subplot(ceil(sub),2,x); 
    y = (1:size(images{x},1));
    plot(y,deltF)
    title(imagefiles(x).name);
    xlabel('Time in Pixels'); % still in px, will convert to ms later
    ylabel('(F-Fo)/Fo');  
    
    clear restings
    clear deltF
end


%% Plot ca2+ flux for first peak of each .tif file
% Peaks are NOT normalized to time = 0 but based on time of actual peak

for x=1:length(images)
    
    % import values
    deltF = deltaF{x};
    TimeMins = TimeMin{x};
    TimeMaxs = TimeMax{x};
    IntMins = IntMin{x};
    IntMaxs = IntMax{x};
    deltF = deltaF{x};
    restings = resting{x};
    
    % loop to get first peak or "n" peak
    for num_peak=1:length(TimeMins)-2;
        
    % first "if" statement is if 'findpeaks()' or 'islocalmin()' improperly found
    % 2 datapoints back to back
        if TimeMins(1,num_peak+1) < TimeMaxs(1,num_peak) 
            valo = TimeMins(1,num_peak+1); % first time minima (i.e. before peak)
            valt = TimeMins(1,num_peak+2); % second time minima (i.e. after peak)
            val_max = TimeMaxs(1,num_peak); % maxima time (i.e. peak)
            deltFo = deltF(1,TimeMins(1,num_peak+1)); % F/F0 before peak
            deltFmax = deltF(1,TimeMaxs(1,num_peak)); % F/F0 at peak

        % second "if" statement is if the first peak is cut off, so the second
        % peak will be selected
        elseif (IntMins(num_peak,1) < (restings - 75)) | (IntMins(num_peak,1) > (restings + 75))
            valo = TimeMins(1,num_peak+1);
            valt = TimeMins(1,num_peak+2);
            val_max = TimeMaxs(1,num_peak+1);
            deltFo = deltF(1,TimeMins(1,num_peak+1)); 
            deltFmax = deltF(1,TimeMaxs(1,num_peak+1));

        % otherwise, peak should be normal
        else
            valo = TimeMins(1,num_peak);
            valt = TimeMins(1,num_peak+1);
            val_max = TimeMaxs(1,num_peak);
            deltFo = deltF(1,TimeMins(1,num_peak)); 
            deltFmax = deltF(1,TimeMaxs(1,num_peak));
        end 
        
        % save values
        valone{num_peak} = valo;
        valtwo{num_peak} = valt;
        val__max{num_peak} = val_max;
        raw_tau_{num_peak} = val_max - valo; % time to peak
        delta_F1{num_peak} = deltFo;
        delta_Fmax{num_peak} = deltFmax;
        delta_F_diff{num_peak} = deltFmax - deltFo;
        
    end
    
    val1{x} = valone;
    val2{x} = valtwo;
    valmax{x} = val__max;
    raw_tau{x} = raw_tau_; % time to peak
    deltaF1{x} = delta_F1;
    deltaFmax{x} = delta_Fmax;
    deltaF_diff{x} = delta_F_diff;
    
%     valone{1} = valone{1} - 50; % add more data to plot
   
    % plot with deltaF values
    sub = length(images)/4; % to make a grid of plots in 4 columns
    subplot(ceil(sub),4,x); % to specify number of rows
    y = valone{1}:valtwo{1};
    
    plot(y*conv,deltF(valone{1}:valtwo{1})) % conversion 'conv' set before this loop
    title(imagefiles(x).name);
    xlabel('Time in Milliseconds');
    ylabel('(F-Fo)/Fo');
    ylim([-0.4 (max(deltF)*1.2)]);
    
    clear valo
    clear valt
    clear val_max
    clear deltFo
    clear deltFmax
    clear valone
    clear valtwo
    clear val__max
    clear raw_tau_
    clear delta_F1
    clear delta_Fmax
    clear delta_F_diff
    clear deltF
    clear TimeMins
    clear TimeMaxs
    clear IntMins
    clear IntMaxs
    clear restings
end

%% Plot ca2+ flux for first peak of each .tif file
% Peaks ARE normalized to time = 0 and saved individually into 'pathout'

for x=1:length(images)
    
    % import values
    deltF = deltaF{x};
    first = val1{x};
    valo = first{1}; % change val for which peak
    second = val2{x};
    valt = second{1}; % change val for which peak
    
%     valo = valo - 50; % add more to plot
    
    % individual plot
    filename = imagefiles(x).name; % extract file name
    filename_split = strsplit(filename,'.'); % remove "."
    final = strcat(filename_split{2},filename_split{1}); % move number to end 
    xx = sprintf('%s.png', final); % final format as png
    
    subplot(1,1,1);
    y = (valo-(valo)):(valt-(valo)); % to start at 0
    
    % for better peaks
    z = smooth(deltF(valo:valt), 0.15, 'loess');
    plot(y*conv,z);
    title(imagefiles(x).name);
    xlabel('Time in Milliseconds');
    ylabel('(F-Fo)/Fo');
    ylim([-0.2 (max(deltF)*1.2)]);
    xlim([0 ((valt-(valo))*conv)])
    
%     saveas(gcf, xx) % save each plot
end

%% Calculate tau value
% tau is time at which F/Fo is 63% to peak

percent = 0.63; % you set yourself

for x=1:length(images)
    
    % import in
    deltF_diff = deltaF_diff{x};
    deltFo = deltaF1{x};
    deltFmax = deltaFmax{x};
    valo = val1{x};
    valt = val2{x};
    val_max = valmax{x};
    deltF = deltaF{x};
    
    % get all tau
    for num_peak=1:length(deltF_diff)-1;
    
        % F tau, flux intensity when tau occurs
        f_tau = (deltF_diff{num_peak}*percent) + deltFo{num_peak}; % F tau val

        % set range for tau
        range = 0.01;

        % time tau
        if valo{1} > val_max{1}
            for j=valo{num_peak}:val_max{num_peak+1}
                if deltF(j) < f_tau + range && deltF(j) > f_tau - range % finding approximation
                    tau = j; % not normalized to time = 0
                    norm_tau = j - valo{num_peak}; % normalized to start of first peak at time = 0
                    lambda = 1/tau;
                    norm_lambda = 1/norm_tau;
                    break
                end
            end
            for j=valo{num_peak}:val_max{num_peak+1}
                if exist('tau', 'var') == 0
                    tau = (val_max{num_peak}-valo{num_peak})*percent + valo{num_peak};
                    norm_tau = tau - valo{num_peak}; % normalized to start of first peak at time = 0
                    lambda = 1/tau;
                    norm_lambda = 1/norm_tau;
                    break
                end
            end
        elseif valo{1} < val_max{1}
            for j=valo{num_peak}:val_max{num_peak}
                if deltF(j) < f_tau + range && deltF(j) > f_tau - range % finding approximation
                    tau = j; % not normalized to time = 0
                    norm_tau = j - valo{num_peak}; % normalized to start of first peak at time = 0
                    lambda = 1/tau;
                    norm_lambda = 1/norm_tau;
                end
            end
            for j=valo{num_peak}:val_max{num_peak+1}
                if exist('tau', 'var') == 0
                    tau = (val_max{num_peak}-valo{num_peak})*percent + valo{num_peak};
                    norm_tau = tau - valo{num_peak}; % normalized to start of first peak at time = 0
                    lambda = 1/tau;
                    norm_lambda = 1/norm_tau;
                end
            end
        end
        
        
        f_taus{num_peak} = f_tau;
        taus{num_peak} = tau;
        taus_ms{num_peak} = tau*conv; % convert to ms
        taus_ms_norm{num_peak} = norm_tau*conv;
        lambdas{num_peak} = lambda;
        lambdas_ms{num_peak} = lambda*conv;
        lambdas_ms_norm{num_peak} = norm_lambda*conv;
        
    end
    
    un = unique(cell2mat(taus_ms_norm));
    mod_taus = isoutlier(un, 'mean',  'ThresholdFactor', 1);
    
    % save values
    allftau{x} = f_taus;
    alltau{x} = taus;
    alltau_ms{x} = taus_ms; % convert to ms
    alltau_ms_norm{x} = un(~mod_taus);
    alllambda{x} = lambdas;
    alllambda_ms{x} = lambdas_ms;
    alllambda_ms_norm{x} = lambdas_ms_norm;
    
    % val you choose to plot
    choice = 1;
    
%     valo{choice} = valo{choice} - 50; % add more to plot
    
    % plot with tau values to make sure they are correct
    sub = length(images)/4; % to make a grid of plots in 4 columns
    subplot(ceil(sub),4,x); % to specify number of rows
    
    hold on
    y = valo{choice}:valt{choice};
    plot(y,deltF(valo{choice}:valt{choice}))
    plot(alltau{x}{choice}, allftau{x}{choice}, '*') % add * where tau is
    
%     % conv
%     plot(y*conv,deltF(valo:valt))
%     plot(tau*conv, f_tau, '*') % add * where tau is
    
    title(imagefiles(x).name);
    xlabel('Time in Milliseconds');
    ylabel('(F-Fo)/Fo');
    ylim([-0.4 (max(deltF)*1.2)]);
    hold off
    
end

%% Find and plot time to peak using all peaks in data

% get peak
for x=1:length(images)

    TimeMins = TimeMin{x};
    TimeMaxs = TimeMax{x};
    
    for n = 1:(length(TimeMaxs)-2)
        if TimeMins(1,1) > TimeMaxs(1,1);
            TimetoPeakss = TimeMaxs(1,n+1) - TimeMins(1,n);
            CycleLength(n) = (TimeMins(1,n+1) - TimeMins(1,n))*conv;
            TimetoPeaks(n) = TimetoPeakss*conv;
        elseif TimeMins(1,1) < TimeMaxs(1,1)
            TimetoPeakss = TimeMaxs(1,n) - TimeMins(1,n);
            CycleLength(n) = (TimeMins(1,n+1) - TimeMins(1,n))*conv;
            TimetoPeaks(n) = TimetoPeakss*conv;
        end
    
        % remove outliers
        Mod_TimetoPeaks = isoutlier(TimetoPeaks, 'mean',  'ThresholdFactor', 1.75);
        TimetoPeak{x} = TimetoPeaks(~Mod_TimetoPeaks);
        
        % remove outliers
        new_cycle = isoutlier(CycleLength, 'mean', 'ThresholdFactor', 1);
        
        Cycle_Length{x} = CycleLength(~new_cycle);
        val_std = std(TimetoPeaks(~Mod_TimetoPeaks));
        val_std_peak{x} = val_std;
        val_error = val_std/length(TimetoPeaks(~Mod_TimetoPeaks));
        val_error_peak{x} = val_error;
        
    end
    
    AvgTimetoPeak = mean(TimetoPeaks(~Mod_TimetoPeaks));
    AvgTimetoPeaks{x} = AvgTimetoPeak;
    AllTimetoPeaks{x} = TimetoPeaks(~Mod_TimetoPeaks);
    AvgCycleLength = mean(CycleLength(~new_cycle));
    AvgCycleLengths{x} = AvgCycleLength;
    AllCycleLength{x} = CycleLength(~new_cycle);
    
    clear TimeMins
    clear TimeMaxs
    clear new_cycle
    clear Mod_TimetoPeaks
    clear TimetoPeakss
    clear CycleLength
    clear TimetoPeaks
end

%% make table with variables and export

% rename files so matlab likes the format
%%% CHANGE based on your file names %%%
% Matlab likes letter then number and will not accept the opposite
for x=1:length(images)
    filename = imagefiles(x).name;
    filename_split = strsplit(filename,'.');
%     filename_split2 = strsplit(filename_split{1},'_');
    filenames{x} = strcat(filename_split{1},filename_split{2});
%     filenames{x} = strcat(filename_split{1});
end

%% Gets all of the values into an excel sheet

A = nan(length(images)*6, 30);
variable_labels = [];
for i=1:length(alltau_ms_norm)
    horiz = horzcat(alltau_ms_norm{i});
    A(i, 1:length(horiz)) = horiz;
    name = 'tau(ms)';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(AllTimetoPeaks)
    horiz = horzcat(AllTimetoPeaks{i});
    new = i+length(alltau_ms_norm);
    A(new, 1:length(horiz)) = horiz;
    name = 'timetopeak(ms)';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(deltaF_diff)
    horiz = horzcat(cell2mat(deltaF_diff{i}));
    new = i+length(alltau_ms_norm)+length(AllTimetoPeaks);
    A(new, 1:length(horiz)) = horiz;
    name = 'deltaF/F0';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(Cycle_Length)
    horiz = horzcat(Cycle_Length{i});
    new = i+length(alltau_ms_norm)+length(AllTimetoPeaks)+length(deltaF_diff);
    A(new, 1:length(horiz)) = horiz;
    name = 'cyclelength(ms)';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(Widths)
    horiz = horzcat(Widths{i});
    new = i+length(alltau_ms_norm)+length(AllTimetoPeaks)+length(deltaF_diff)+length(Cycle_Length);
    A(new, 1:length(horiz)) = horiz;
    name = 'widths(ms)';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(Proms)
    horiz = horzcat(Proms{i});
    new = i+length(alltau_ms_norm)+length(AllTimetoPeaks)+length(deltaF_diff)+length(Cycle_Length)+length(Widths);
    A(new, 1:length(horiz)) = horiz;
    name = 'prominance';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end

final_tab = array2table(A);
new_table = [ table(transpose(variable_labels), 'VariableNames',{'Labels'})  final_tab];
writetable(new_table, 'T14_D50_9_75_10k_1000ms.csv');

