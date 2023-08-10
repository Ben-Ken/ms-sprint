function [] = Broadband_periodic_sensitivities(filename, cf_mode, parameter, type)
%BROADBAND_PERIODIC_SENSITIVITIES Plots the peak detection rate of a given
%periodic parameter as a function of frequency
%   Arguments:: Filename: file to import, cf_mode: 'bw' or 'rgb', parameter: 'cf',
%   'amp' or 'sd', type: 'PP' or 'MS'
    
    %Load file in
    load(['Results/' filename '10000.mat'])

    %Create a time point matrix
    ts(:) = linspace(1.5, 58.5, 115);
    
    %Extract relevant data
    c = 1;
    for s = 1:length(vstmp)
        for i = 1:size(vstmp(s).pk_cf,1)
        c = c+1;
        err_timefreq(c,1) = ostmp(s).pk_times(i,2); % duration, in s
        err_timefreq(c,2) = ostmp(s).pk_cf(i); % centre frequency
        err_timefreq(c,3) = ostmp(s).pk_sd(i); % sd
        err_timefreq(c,4) = ostmp(s).pk_amp(i); % amplitude
        err_timefreq(c,5) = sum(~isnan(rstmp(s).pk_cf(i,:)))./sum(~isnan(vstmp(s).pk_cf(i,:))); % peak detection rate
        end
    end
    
    %Set swatch
    switch type
        case 'PP'
            swatch = [0,0.447,0.698];
        case'MS'
            swatch = [0.835,0.369,0];
    end

    switch parameter
        case 'cf' %peak detection rate by center frequency
            switch cf_mode  %alternatives are 'bw' (black and white) or 'rgb' (color)
                case 'bw'
                    scatter(err_timefreq(:,2),err_timefreq(:,5),6,swatch,'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
                case 'rgb' %splits up data into frequnecy bands, each with their own color, splitting at 8Hz, 13Hz, 18Hz
                    pars = [0.85 4];
                    q = 0;
                    swatch = [108 156 167;236 189 78;198 97 43;162 55 27;91 61 28;86 34 98]./255;
                    hold on
                    scatter((err_timefreq((err_timefreq(:,2)<8),2)),err_timefreq((err_timefreq(:,2)<8),5),6,swatch(1,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
                    scatter((err_timefreq((8<=err_timefreq(:,2) & err_timefreq(:,2)<13),2)),err_timefreq((8<=err_timefreq(:,2) & err_timefreq(:,2)<13),5),6,swatch(2,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
                    scatter((err_timefreq((13<=err_timefreq(:,2) & err_timefreq(:,2)<18),2)),err_timefreq((13<=err_timefreq(:,2) & err_timefreq(:,2)<18),5),6,swatch(3,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
                    scatter((err_timefreq((18<=err_timefreq(:,2)),2)),err_timefreq(18<=(err_timefreq(:,2)),5),6,swatch(4,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
                    hold off
            end
            xlabel('Center Frequency (Hz)')
            ylabel('p(peak)')
        case 'amp' %peak detection rate by peak amplitude
            scatter(err_timefreq(:,4),err_timefreq(:,5),6,swatch,'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
            xlabel('Amplitude (a.u.)')
            ylabel('p(peak)')
            xlim([0.6 1.6])
    
        case 'sd' %peak detection rate by standard deviation
            scatter(err_timefreq(:,3),err_timefreq(:,5),6,swatch,'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
            xlabel('Standard Deviation (Hz)')
            ylabel('p(peak)')
            xlim([1 2])
    end
end
