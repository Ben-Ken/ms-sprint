function [] = Broadband_MAE()
%BROADBAND_MAE Plots Mean Absolute Error as a function of frequency across
%models

    hold on
    files = {'conservative' 'MSnatural2' 'MSconservative'};
    name = {'Conservative' 'MS' 'MS-conservative'};
    col = [0.337,0.706,0.914,0,0.620,0.451,0,0.447,0.698];
    legend('FontSize', 6)

    %for every model
    for i = 1:3
        load(['Results/SPECTRA' files{i} '10000.mat']) % load file
        freq_error = nan(length(spectra),40); %pre-allocate
        for sim = 1:length(spectra)
            freq_error(sim, :) = spectra(sim).frequency_wise_error; %import freq-wise error
        end
        freq_error = freq_error'; %vectorize
        
        %taking the average
        freq_error = mean(freq_error, 2);
        plot(1:40,freq_error, 'DisplayName', name{i},'Color',[col(((i-1)*3+1):i*3)], 'LineWidth', 1.5)
        ylim([0.1 0.2])
    end
end
