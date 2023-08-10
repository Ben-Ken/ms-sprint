function [] = Aperiodic_error_bypk()
%APERIODIC_ERROR_BYPK Plots boxchart of MAE of exponent for every number of simulated peaks by model 
    
    %Chosen models and colors
    files = {'conservative' 'MSnatural2' 'MSconservative'};
    name = {'Conservative' 'MS' 'MS-conservative'};
    nchan = 10000;
    col = [0.337,0.706,0.914,0,0.620,0.451,0,0.447,0.698]; %set swatch
    hold on
    for t = 1:3 %for every model
        clear per_peak_mat
        load(['Results/' files{t} '10000.mat']) %import file

        %Calculate absolute error, reshape into chans*timepoints
        expE = abs([rstmp.ap_exp] - [vstmp.ap_exp]);
        expE = reshape(expE, [115,nchan]).';

        %number of peaks per channel
        expE(:,116) = [rstmp(:).n_peaks];

        %for every number of peaks in a channel, extract all of the errors
        per_peak_mat = nan(5,1150000); %preallocate
        for n = 0:max(expE(:,116)) %for every number of peaks
            clear per_peak_tmp %clear previous
            rows = size([expE((expE(:,116) == n), :)],1); %number of rows with n peaks
            per_peak_tmp(1:rows,1:115) = [expE((expE(:,116) == n), 1:115)]; %create matrix of errors w/ n peaks
            per_peak_mat(n+1,1:length([per_peak_tmp(:)])) = per_peak_tmp(:); %vectorize into nth row of matrix
        end
        per_peak_mat(per_peak_mat == 0) = NaN; %filter out products of varying amounts of data
        per_peak_mat = per_peak_mat';

        %Create boxchart in relevant column grouping by peak #, ordered and
        %colored by model
        column_mat = nan(1150000,15);
        column_mat(:,t) = per_peak_mat(:,1);
        column_mat(:,3+t) = per_peak_mat(:,2);
        column_mat(:,6+t) = per_peak_mat(:,3);
        column_mat(:,9+t) = per_peak_mat(:,4);
        column_mat(:,12+t) = per_peak_mat(:,5);
        boxchart(column_mat, 'BoxFaceColor', col(1,(3*t-2:3*t)), 'MarkerColor','#FFFFFF', 'MarkerSize',0.01);
        ylim([0 0.5])
    end

    %Clean up figure
    xticklabels(["" "0 peaks" "" "" "1 peak" "" "" "2 peaks" "" "" "3 peaks" "" "" "4 peaks" ""])
    xtickangle(0)
end


