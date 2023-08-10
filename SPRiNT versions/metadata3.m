%Calculation of PPV and sensitivity by taking the sum and dividing by the
%total. This gives data *by peak*
function [stats, stats_mats] = metadata3(mAE_cell,vstmp,rstmp)

stats_mats = nan(length(rstmp),5);
% TP, # sim peaks, # modelled peaks, by sim sens, by sim PPV

for s = 1:length(rstmp)
    sum_exp = 0;
    if ~(rstmp(s).n_peaks == 0)
        for pk = 1:length(rstmp(s).exp_peaks)
            sum_exp = sum_exp + length(rstmp(s).exp_peaks(pk).peak);
        end
        TP = sum_exp;
        stats_mats(s,1) = TP;
        FP = length(rstmp(s).resid_peaks);
        stats_mats(s,2) = sum(vstmp(s).pk_exist, 'all');     %equivalent to TP + FN
        stats_mats(s,3) = TP + FP;
        stats_mats(s,5) = TP/(TP+FP); %PPV
    else
        stats_mats(s,3) = length(rstmp(s).resid_peaks);
    end 
    stats_mats(s,4) = TP/(stats_mats(s,2)); %sens
end

stats.sens = sum(stats_mats(:,1),'omitnan')/sum(stats_mats(:,2),'omitnan');
stats.PPV = sum(stats_mats(:,1),'omitnan')/sum(stats_mats(:,3),'omitnan');

CImat = bootstrap3(stats_mats, 5000);
stats.sens_CI = CImat(1,:);
stats.PPV_CI = CImat(2,:);

%F1 score :)
stats.F1 = 2/((stats.sens)^-1+(stats.PPV)^-1);
stats.F1_CI = [sqrt(((1/(stats.sens)+1/(stats.PPV))^(-2)*(2/(stats.sens)^2)*mean((stats.sens_CI(1))))^2+...
    ((1/(stats.PPV)+1/(stats.sens))^(-2)*(2/(stats.PPV)^2)*mean((stats.PPV_CI(1))))^2),...
               sqrt(((1/(stats.sens)+1/(stats.PPV))^(-2)*(2/(stats.sens)^2)*mean((stats.sens_CI(2))))^2+...
    ((1/(stats.PPV)+1/(stats.sens))^(-2)*(2/(stats.PPV)^2)*mean((stats.PPV_CI(2))))^2)];

% mMAE for aperiodic paramters
stats.apexp_mMAE = mean(10.^(mAE_cell{:,1}));
stats.apoff_mMAE = mean(10.^(mAE_cell{:,2}));

%mMAE for periodic parameters
stats.cf_mMAE = mean(10.^(mAE_cell{:,3}));
stats.amp_mMAE = mean(10.^(mAE_cell{:,4}));
stats.sd_mMAE = mean(10.^(mAE_cell{:,5}));



