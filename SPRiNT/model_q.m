function [stats,stats_mat, mAE_cell, rstmp, vstmp, ostmp] = model_q(s_data, chans)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Load data
os = open("generic_simulation_parameters.mat");
os = os.generic_simulation_parameters(1:chans);

SPRiNT = s_data;
work = whos('SPRiNT'); % Data is saved as SPRiNT01, SPRiNT02, etc.
if numel(work) > 1
    sd = eval(work(1).name).SPRiNT;
    for t = 2:numel(work)
        sdtmp = eval(work(t).name).SPRiNT;
        sd.channel(end+1:end+length(sdtmp.channel)) = sdtmp.channel;
    end
else
    sd = eval(work(1).name).SPRiNT;
end

%
ts = [sd.channel(1).data.time];
len = 60; % length of recording, in seconds
fs = 200;
tvec = 1/fs:1/fs:60;
tuk = tukeywin(1000,0.4);

sdtmp = sd; % SPRiNT data structures
ostmp = os; % Simulation parameters
vstmp = struct(); % Expected output statistics (under perfect conditions)
rstmp = struct(); % SPRiNT output statistics

% Generate expected peak statistics
for sim = 1:length(sdtmp.channel)
    vstmp(sim).ap_exp = zeros(fs.*len,1);
    vstmp(sim).ap_off = zeros(fs.*len,1);
    
    vstmp(sim).ap_exp(1:round(ostmp(sim).ap_dyn_prop(1).*len.*fs)) = ostmp(sim).ap_init(1);
    vstmp(sim).ap_exp(round(ostmp(sim).ap_dyn_prop(1).*len.*fs)+1:round(ostmp(sim).ap_dyn_prop(2).*len.*fs)) = ...
        linspace(ostmp(sim).ap_init(1),ostmp(sim).ap_init(1)+ostmp(sim).ap_chng(1),round(ostmp(sim).ap_dyn_prop(2).*len.*fs)-round(ostmp(sim).ap_dyn_prop(1).*len.*fs));
    vstmp(sim).ap_exp(round(ostmp(sim).ap_dyn_prop(2).*len.*fs)+1:end) = ostmp(sim).ap_init(1)+ostmp(sim).ap_chng(1);
    vstmp(sim).ap_exp = spline(tvec,vstmp(sim).ap_exp,ts);
    
    vstmp(sim).ap_off(1:round(ostmp(sim).ap_dyn_prop(1).*len.*fs)) = ostmp(sim).ap_init(2);
    vstmp(sim).ap_off(round(ostmp(sim).ap_dyn_prop(1).*len.*fs)+1:round(ostmp(sim).ap_dyn_prop(2).*len.*fs)) = ...
        linspace(ostmp(sim).ap_init(2),ostmp(sim).ap_init(2)+ostmp(sim).ap_chng(2),round(ostmp(sim).ap_dyn_prop(2).*len.*fs)-round(ostmp(sim).ap_dyn_prop(1).*len.*fs));
    vstmp(sim).ap_off(round(ostmp(sim).ap_dyn_prop(2).*len.*fs)+1:end) = ostmp(sim).ap_init(2)+ostmp(sim).ap_chng(2);
    vstmp(sim).ap_off = spline(tvec,vstmp(sim).ap_off,ts);
    
    vstmp(sim).pk_cf = nan(size(ostmp(sim).pk_cf,1),length(ts));
    vstmp(sim).pk_amp = nan(size(ostmp(sim).pk_cf,1),length(ts));
    vstmp(sim).pk_sd = nan(size(ostmp(sim).pk_cf,1),length(ts));
    vstmp(sim).pk_exist = zeros(size(ostmp(sim).pk_cf,1),length(ts));
    
    for pk = 1:size(ostmp(sim).pk_cf,1)
        vstmp(sim).pk_exist(pk,:) = (ostmp(sim).pk_times(pk,1) <= ts) & (ts <= ostmp(sim).pk_times(pk,1)+ostmp(sim).pk_times(pk,2));
        vstmp(sim).pk_cf(pk,logical(vstmp(sim).pk_exist(pk,:))) = ostmp(sim).pk_cf(pk);
        vstmp(sim).pk_amp(pk,logical(vstmp(sim).pk_exist(pk,:))) = spline(1:1000,tuk,linspace(1,1000,sum(vstmp(sim).pk_exist(pk,:)))).*ostmp(sim).pk_amp(pk);
        vstmp(sim).pk_sd(pk,logical(vstmp(sim).pk_exist(pk,:))) = ostmp(sim).pk_sd(pk);
    end
    
    rstmp(sim).n_peaks = size(ostmp(sim).pk_cf,1);
    
end

% Extract peak statistics
sdtmp = sd;
for sim = 1:length(sdtmp.channel)
    npeak_vec = zeros(1,length(ts));
    npeak_exp = sum(~isnan(vstmp(sim).pk_cf),1);
    if isempty(npeak_exp)
        npeak_exp = zeros(1,length(ts));
    end
    peak_props = zeros(5,7,6);
    for i = 1:length(ts)
        npeak_vec(i) = sum(round([sdtmp.channel(sim).peaks.time],1)==round(ts(i),1));
    end
    for i = 1:length(ts)
        pks_tmp = [sdtmp.channel(sim).peaks(round([sdtmp.channel(sim).peaks.time],1)==round(ts(i),1)).center_frequency];
        for j = 1:length(pks_tmp)
            if pks_tmp(j) < 8
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,1) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,1)+1;
            elseif 8 <= pks_tmp(j) && pks_tmp(j) < 13;
                %gives me an error sometimes because of the 's'
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,2) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,2)+1;
            elseif 13 <= pks_tmp(j) && pks_tmp(j) < 18
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,3) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,3)+1;
            elseif 18 <= pks_tmp(j) && pks_tmp(j) < 23
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,4) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,4)+1;
            elseif 23 <= pks_tmp(j) && pks_tmp(j) < 28
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,5) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,5)+1;
            elseif 28 <= pks_tmp(j) && pks_tmp(j) < 33
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,6) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,6)+1;
            else
            end
        end
    end
    peaks = sdtmp.channel(sim).peaks;
    keep = zeros(1,length(peaks));
    %%%%%% This represents a departure from how we have been doing peak
    %%%%%% detection
    zn = [ostmp(sim).pk_cf-2.5.*ostmp(sim).pk_sd ostmp(sim).pk_cf+2.5.*ostmp(sim).pk_sd]; %because the ms-specparam paper used 2.5 sd, 
    if size(ostmp(sim).pk_cf,1)                                                           % I switched the 2 -> 2.5 here, I don't know 
        for pk = 1:size(ostmp(sim).pk_cf,1)                                               %if there were other changes to make elsewhere
            keep = keep | ((ostmp(sim).pk_times(pk,1) <= [peaks.time]) & ([peaks.time] <= ostmp(sim).pk_times(pk,1)+ostmp(sim).pk_times(pk,2)) & (zn(pk,1) <= [peaks.center_frequency]) & (zn(pk,2) >= [peaks.center_frequency]));
        end
    else
        keep = zeros(1,length(peaks));
    end
    discard = ~keep;
    rstmp(sim).resid_peaks = peaks(discard);
    sdtmp.channel(sim).peaks(discard) = [];
    peaks(discard) = [];    
    if size(ostmp(sim).pk_cf,1)
        [tmp, order] = sort(ostmp(sim).pk_amp,'descend');
        rstmp(sim).exp_peaks = struct();
        for ind = 1:length(order)
            pk = order(ind);
            expi = find((ostmp(sim).pk_times(pk,1) <= [peaks.time]) &...
                ([peaks.time] <=...
                ostmp(sim).pk_times(pk,1)+ostmp(sim).pk_times(pk,2)) &...
                (zn(pk,1) <= [peaks.center_frequency]) & (zn(pk,2) >=...
                [peaks.center_frequency])); % filter by cf and time range
%             expi = find((zn(pk,1) <= [peaks.center_frequency]) & (zn(pk,2) >= [peaks.center_frequency])); % filter by cf range only
            exp_peaks = peaks(expi);
            if length(exp_peaks) > length(unique([exp_peaks.time]))
                for i = unique([exp_peaks.time])
                    pks = find([exp_peaks.time] == i);
                    if length(pks) >1 
                        [tmp, rmv] = max([exp_peaks([exp_peaks.time] == i).amplitude]);
                        pks(rmv) = [];
                        exp_peaks(pks) = [];
                        expi(pks) = [];
                    end   
                end
            end
            rstmp(sim).exp_peaks(pk).peak = peaks(expi);
            peaks(expi) = [];
        end
        if ~isempty(peaks)
            rstmp(sim).resid_peaks(length(rstmp(sim).resid_peaks)+1:length(rstmp(sim).resid_peaks)+length(peaks)) = peaks;
        end
    end

    rstmp(sim).ap_exp = [sdtmp.channel(sim).aperiodics.exponent];
    rstmp(sim).ap_off = [sdtmp.channel(sim).aperiodics.offset];
    rstmp(sim).pk_cf = nan(size(vstmp(sim).pk_cf));
    rstmp(sim).pk_amp = nan(size(vstmp(sim).pk_amp));
    rstmp(sim).pk_sd = nan(size(vstmp(sim).pk_sd));
    rstmp(sim).pk_det = nan(size(vstmp(sim).pk_cf));
    
    for pk = 1:size(ostmp(sim).pk_cf,1)
        [~, iMat, iPk] = intersect(ts,[rstmp(sim).exp_peaks(pk).peak.time]);
        rstmp(sim).pk_cf(pk,iMat) = [rstmp(sim).exp_peaks(pk).peak(iPk).center_frequency];
        rstmp(sim).pk_amp(pk,iMat) = [rstmp(sim).exp_peaks(pk).peak(iPk).amplitude];
        rstmp(sim).pk_sd(pk,iMat) = [rstmp(sim).exp_peaks(pk).peak(iPk).st_dev];
        rstmp(sim).pk_det(pk,~isnan(vstmp(sim).pk_cf(pk,:))) = 0;
        rstmp(sim).pk_det(pk,iMat) = 1;
    end
    
    rstmp(sim).expmod_peaks = [npeak_exp;npeak_vec];
    rstmp(sim).peak_props = peak_props;
    
end
%%
% Calculate sensitivity and specificity
sens_mat = zeros(length(rstmp),2);
spec_mat = zeros(length(rstmp),1);
for i = 1:length(rstmp)
    sens_mat(i,1) = sum(vstmp(i).pk_exist(:)); % every peak that exists in the given channel
    sens_mat(i,2) = sum(rstmp(i).pk_det(logical(vstmp(i).pk_exist)));
    spec_mat(i) = length(rstmp(i).resid_peaks);
end
sens = sum(sens_mat(:,2))./sum(sens_mat(:,1));
spec = 1-sum(spec_mat)./length([sdtmp.channel(1).data.time])./length(sdtmp.channel);

% Extract a variety of peak descriptives and errors
err_timefreq = nan(40000,10);
c = 0;
for s = 1:length(sd.channel)
    for i = 1:size(vstmp(s).pk_cf,1)
        c = c+1;
        err_timefreq(c,1) = ostmp(s).pk_times(i,2); % duration, in s
        err_timefreq(c,2) = ostmp(s).pk_cf(i); % centre frequency
        err_timefreq(c,3) = mean(abs(rstmp(s).pk_cf(i,:)-vstmp(s).pk_cf(i,:)),'omitnan'); % cf MAE
        err_timefreq(c,4) = mean(abs(rstmp(s).pk_amp(i,:)-vstmp(s).pk_amp(i,:)),'omitnan'); % amp MAE
        err_timefreq(c,5) = mean(abs(rstmp(s).pk_sd(i,:)-vstmp(s).pk_sd(i,:)),'omitnan'); % sd MAE
        err_timefreq(c,6) = ostmp(s).pk_sd(i); % sd
        err_timefreq(c,7) = ostmp(s).pk_times(i,2).*ostmp(s).pk_cf(i); % ncycles
        err_timefreq(c,8) = sum(~isnan(rstmp(s).pk_cf(i,~isnan(vstmp(s).pk_cf(i,:)))))./sum(~isnan(vstmp(s).pk_cf(i,:))); % peak detection rate
        err_timefreq(c,9) = mean(vstmp(s).ap_exp(~isnan(vstmp(s).pk_cf(i,:)))); % mean aperiodic slope
        err_timefreq(c,10) = ostmp(s).pk_amp(i); % amplitude
    end
end
err_timefreq(c+1:end,:) = [];

% Calculate error by number of peaks
err_bypk = zeros(3,5);
err_bypk(1,:) = 0:4;
err_bypkmat = zeros(1150000,5);
c = 0;  
for s = 1:length(sd.channel)

    for i = 1:115
        c = c+1;
        if isempty(vstmp(s).pk_cf)
            err_bypkmat(c,1) = 0;
        else
            err_bypkmat(c,1) = sum(~isnan(vstmp(s).pk_cf(:,i)));
        end
        err_bypkmat(c,2) = sqrt(40.*sdtmp.channel(s).stats(i).MSE)./40; % MAE
    end
end
err_bypkmat(c+1:end,:) =[];
err_bypkcell = {err_bypkmat(err_bypkmat(:,1) == 0,2),err_bypkmat(err_bypkmat(:,1) == 1,2),...
    err_bypkmat(err_bypkmat(:,1) == 2,2),err_bypkmat(err_bypkmat(:,1) == 3,2),err_bypkmat(err_bypkmat(:,1) == 4,2)};
% calculate error for each peak parameter
err_peaks = zeros(400000,3);
err_aper = zeros(length(sd.channel).*length(sd.channel(1)),2);
ctr = 1; 
for s = 1:length(sd.channel)
    for i = 1:size(vstmp(s).pk_cf,1)
        a = abs(rstmp(s).pk_cf(i,:)-vstmp(s).pk_cf(i,:));
        b = abs(rstmp(s).pk_amp(i,:)-vstmp(s).pk_amp(i,:));
        c = abs(rstmp(s).pk_sd(i,:)-vstmp(s).pk_sd(i,:));
        idx = ~isnan(a);
        err_peaks(ctr:ctr+sum(idx)-1,1) = a(idx);
        err_peaks(ctr:ctr+sum(idx)-1,2) = b(idx);
        err_peaks(ctr:ctr+sum(idx)-1,3) = c(idx);
        ctr = ctr+sum(idx);
    end
    err_aper((s-1)*length(ts)+1:s*length(ts),1) = abs(rstmp(s).ap_off-vstmp(s).ap_off);
    err_aper((s-1)*length(ts)+1:s*length(ts),2) = abs(rstmp(s).ap_exp-vstmp(s).ap_exp);
end
err_peaks(ctr:end,:) = [];
err_peaks = log10(err_peaks);
err_aper = log10(err_aper);
mAE_cell = {err_aper(:,2),err_aper(:,1),err_peaks(:,1),err_peaks(:,2),err_peaks(:,3)}; % Exponent, offset, cf, amp,sd

npeak_mat = zeros(5,7);
for s = 1:length(sdtmp.channel)
    for i = 1:length([rstmp(s).expmod_peaks])
        npeak_mat(rstmp(s).expmod_peaks(1,i)+1,rstmp(s).expmod_peaks(2,i)+1) =...
            npeak_mat(rstmp(s).expmod_peaks(1,i)+1,rstmp(s).expmod_peaks(2,i)+1)+1;
    end
end
pp = zeros(5,7,6);
for i = 1:length(sdtmp.channel)
    pp = pp + rstmp(i).peak_props;
end
% Focus on first three freq ranges
pp(:,:,4:6) = [];

%TEMPORARY FUNCTION FOR PRELIMENARY RESULTS
[stats, stats_mat] = metadata3(mAE_cell,vstmp,rstmp);

end

