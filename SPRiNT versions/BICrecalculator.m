%the BIC recalculator (load SPRiNT and fs into workspace and let it work
%its magic!
%=========================================================================
function [SPRiNT] = BIC_recalculator(SPRiNT, fs)
    switch SPRiNT.options.peak_type 
        case 'gaussian' % gaussian only
            peak_function = @gaussian;
        case 'cauchy'
            peak_function = @cauchy;
    end
    timeRange = SPRiNT.options.maxtime.*SPRiNT.options.winLen.*(1-SPRiNT.options.Ovrlp./100);
    nC = length(SPRiNT.channel);
    for c = 1:nC
        ts = [SPRiNT.channel(c).data.time];
        for t = 1:length(ts)
            peak_fit = zeros(size(SPRiNT.freqs));
            p = [SPRiNT.channel(c).peaks.time] == ts(t);
            SPRiNT.channel(c).data(t).peak_params = [[SPRiNT.channel(c).peaks(p).center_frequency]' [SPRiNT.channel(c).peaks(p).amplitude]' [SPRiNT.channel(c).peaks(p).st_dev]'];
            if isempty(SPRiNT.channel(c).data(t).peak_params)
                SPRiNT.channel(c).data(t).peak_params = [0,0,0];
            end
            peak_pars = SPRiNT.channel(c).data(t).peak_params;
            %fit peak_fit with ever peak at this time point 
            for peak = 1:size(peak_pars,1)
                if ~(SPRiNT.channel(c).data(t).peak_params(peak,:) == 0)
                    peak_fit = peak_fit + peak_function(SPRiNT.freqs,peak_pars(peak,1),...
                        peak_pars(peak,2),peak_pars(peak,3));
                end
            end
            ap_spec = log10(SPRiNT.channel(c).data(t).power_spectrum) - peak_fit;
            ap_pars = simple_ap_fit(SPRiNT.freqs, ap_spec, SPRiNT.options.aperiodic_mode, SPRiNT.channel(c).data(t).aperiodic_params(end));
            ap_fit = gen_aperiodic(SPRiNT.freqs, ap_pars, SPRiNT.options.aperiodic_mode);
            MSE = sum((ap_spec - ap_fit).^2)/length(SPRiNT.freqs);
            rsq_tmp = corrcoef(ap_spec+peak_fit,ap_fit+peak_fit).^2;
            % Return FOOOF results
            ap_pars(2) = abs(ap_pars(2));
            SPRiNT.channel(c).data(t).ap_fit = 10.^(ap_fit);
            SPRiNT.channel(c).data(t).fooofed_spectrum = 10.^(ap_fit+peak_fit);
            SPRiNT.channel(c).data(t).peak_fit = 10.^(peak_fit);
            SPRiNT.aperiodic_models(c,t,:) = SPRiNT.channel(c).data(t).ap_fit;
            SPRiNT.SPRiNT_models(c,t,:) = SPRiNT.channel(c).data(t).fooofed_spectrum;
            SPRiNT.peak_models(c,t,:) = SPRiNT.channel(c).data(t).peak_fit;
            SPRiNT.channel(c).aperiodics(t).offset = ap_pars(1);
            if length(ap_pars)>2 % Legacy FOOOF alters order of parameters
                SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(3);
                SPRiNT.channel(c).aperiodics(t).knee_frequency = ap_pars(2);
            else
                SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(2);
            end
            SPRiNT.channel(c).stats(t).MSE = MSE;
            SPRiNT.channel(c).stats(t).r_squared = rsq_tmp(2);
            SPRiNT.channel(c).stats(t).frequency_wise_error = abs(ap_spec-ap_fit);

            %%%% BIC, BF module (by Ben)
            model_fit = ap_fit + peak_fit;
            loglik = -length(model_fit)/2.*(1+log(SPRiNT.channel(c).stats(t).MSE)+log(2*pi));
            SPRiNT.channel(c).data(t).loglik = loglik;
            if SPRiNT.channel(c).data(t).peak_params(1) == 0
                BIC = length(ap_pars).*log(length(model_fit))-2*loglik;
            else
                BIC = (length(ap_pars) + numel(SPRiNT.channel(c).data(t).peak_params)).*log(length(model_fit))-2*loglik;
            end
            SPRiNT.channel(c).data(t).BIC = BIC;
            BF = exp((BIC-SPRiNT.channel(c).data(t).models(1).BIC)./2);
            %%%%
            end
        end
end
  
    
%SUPPORT FUNCTIONS    
   
%% ===== GENERATE APERIODIC =====
function ap_vals = gen_aperiodic(freqs,aperiodic_params,aperiodic_mode)
%       Generate aperiodic values, from parameter definition.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%       	Frequency vector to create aperiodic component for.
%       aperiodic_params : 1x3 array
%           Parameters that define the aperiodic component.
%       aperiodic_mode : {'fixed', 'knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       ap_vals : 1d array
%           Generated aperiodic values.

    switch aperiodic_mode
        case 'fixed'  % no knee
            ap_vals = expo_nk_function(freqs,aperiodic_params);
        case 'knee'
            ap_vals = expo_function(freqs,aperiodic_params);
        case 'floor'
            ap_vals = expo_fl_function(freqs,aperiodic_params);
    end
end

%% ===== CORE MODELS =====
function ys = gaussian(freqs, mu, hgt, sigma)
%       Gaussian function to use for fitting.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create gaussian fit for.
%       mu, hgt, sigma : doubles
%           Parameters that define gaussian function (centre frequency,
%           height, and standard deviation).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for gaussian function.

    ys = hgt*exp(-(((freqs-mu)./sigma).^2) /2);

end

function ys = cauchy(freqs, ctr, hgt, gam)
%       Cauchy function to use for fitting.
% 
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create cauchy fit for.
%       ctr, hgt, gam : doubles
%           Parameters that define cauchy function (centre frequency,
%           height, and "standard deviation" [gamma]).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for cauchy function.

    ys = hgt./(1+((freqs-ctr)/gam).^2);

end

function ys = expo_function(freqs,params)
%       Exponential function to use for fitting 1/f, with a 'knee' (maximum at low frequencies).
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       params : 1x3 array (offset, knee, exp)
%           Parameters (offset, knee, exp) that define Lorentzian function:
%           y = 10^offset * (1/(knee + x^exp))
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential function.

    ys = params(1) - log10(abs(params(2)) +freqs.^params(3));

end

function ys = expo_nk_function(freqs, params)
%       Exponential function to use for fitting 1/f, without a 'knee'.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       params : 1x2 array (offset, exp)
%           Parameters (offset, exp) that define Lorentzian function:
%           y = 10^offset * (1/(x^exp))
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential (no-knee) function.

    ys = params(1) - log10(freqs.^params(2));

end

function ys = expo_fl_function(freqs, params)

    ys = log10(f.^(params(1)) * 10^(params(2)) + params(3));

end

%% ===== FITTING ALGORITHM =====
function aperiodic_params = simple_ap_fit(freqs, power_spectrum, aperiodic_mode, aperiodic_guess)
%       Fit the aperiodic component of the power spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       aperiodic_mode : {'fixed','knee'}
%           Defines absence or presence of knee in aperiodic component.
%       aperiodic_guess: double
%           SPRiNT specific - feeds previous timepoint aperiodic slope as
%           guess
%
%       Returns
%       -------
%       aperiodic_params : 1xn array
%           Parameter estimates for aperiodic fit.

%       Set guess params for lorentzian aperiodic fit, guess params set at init
    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);

    switch (aperiodic_mode)
        case 'fixed'  % no knee
            guess_vec = [power_spectrum(1), aperiodic_guess];
            aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs, power_spectrum);
        case 'knee'
            guess_vec = [power_spectrum(1),0, aperiodic_guess];
            aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs, power_spectrum);
    end

end

%% ===== ERROR FUNCTIONS =====
function err = error_expo_nk_function(params,xs,ys)
    ym = -log10(xs.^params(2)) + params(1);
    err = sum((ys - ym).^2);
end

function err = error_expo_function(params,xs,ys)
    ym = expo_function(xs,params);
    err = sum((ys - ym).^2);
end

function err = error_model(params, xVals, yVals, peak_type, guess, guess_weight)
    fitted_vals = 0;
    weak = 1E2;
    strong = 1E7;
    for set = 1:size(params,1)
        switch (peak_type)
            case 1 % Gaussian
                fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
            case 2 % Cauchy
                fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
        end
    end
    switch guess_weight
        case 'none'
            err = sum((yVals - fitted_vals).^2);
        case 'weak' % Add small weight to deviations from guess m and amp
            err = sum((yVals - fitted_vals).^2) + ...
                 weak*sum((params(:,1)-guess(:,1)).^2) + ...
                 weak*sum((params(:,2)-guess(:,2)).^2);
        case 'strong' % Add large weight to deviations from guess m and amp
            err = sum((yVals - fitted_vals).^2) + ...
                 strong*sum((params(:,1)-guess(:,1)).^2) + ...
                 strong*sum((params(:,2)-guess(:,2)).^2);
    end
end