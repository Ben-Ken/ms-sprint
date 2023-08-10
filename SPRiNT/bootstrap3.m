function [CImat] = btstp(data, sims)

    % sims = 5000;
    % buckets = 35;
    % data = stats_mats;
    
    for line = 1:size(data,1)
        clear(line) = isnan(data(line,1));
    end
    data(clear,:) = [];
    
    %sens
    t_sens = sum(data(:,1),'omitnan')/sum(data(:,2),'omitnan');
    
    for p = 1:sims
        index = randi(size(data,1), 1, size(data,1));
        %create sens based off of random (index) peaks and simulated's
        sens_b(p) = sum(data(index,1),'omitnan')/sum(data(index,2),'omitnan');
    end
    CI = [prctile(sens_b, 2.5) prctile(sens_b, 97.5)];
    
    minCI = abs(t_sens - CI(1));
    maxCI = abs(t_sens - CI(2));
    CImat(1,:) = [minCI, maxCI];
    
    %PPV
    t_PPV = sum(data(:,1),'omitnan')/sum(data(:,3),'omitnan');
    
    for p = 1:sims
        index = randi(size(data,1), 1, size(data,1));
        %create PPV based off of random (index) peaks and simulated's
        PPV_b(p) = sum(data(index,1),'omitnan')/sum(data(index,3),'omitnan');
    end
    CI = [prctile(PPV_b, 2.5) prctile(PPV_b, 97.5)];
    
    minCI = abs(t_PPV - CI(1));
    maxCI = abs(t_PPV - CI(2));
    CImat(2,:) = [minCI, maxCI];

end

% %plotting
% figure
% subplot(1,2,1)
% histogram(sens_b, buckets);
% xlim([0.5 1])
% 
% 
% subplot(1,2,2)
% 
% hold on 
% bar(1, t_mean);
% errorbar(1, t_mean, minCI, maxCI, '.');
% ylim([0,1])
% 
% hold off
% %end