clear all
close all
clc

rng('shuffle');

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/ModelV1Uncertainty/EstimationTask/Scripts/')

% TIP: Keep all possible 1D variables to row vector
% Note: Time step should be set carefully. It should be pretty small
% relative to the spike rate of single neuron so that the probability of
% firing does not shoot up.
% TODO: add a check to make sure firing rate is not so big compared to the
% time window.

% ----------------------------------
% Params
% ----------------------------------
nNeurons = 200;      % Count of neurons
stimDuration = 1;    % Stimulus duration set to 1 seconds
timeStep = 0.001;    % 0.001s (1ms) - Step size of time bins used for binning stimulus duration 
% ntrials = 100;     % Offloaded to generateStim function. No of stimulus to simulate the model with

neuronsPrefOrientation = zeros(nNeurons, 1);
timeBins               = 0:timeStep:stimDuration; 
stimRespProfile        = 1 + zeros(1, numel(timeBins));

% Neurons tuning parameters
tuningParams.d     = zeros(nNeurons, 1) + 0;   % (fixed) Direction selectivity - set it to zero (no need for neuron to be directional selective).
tuningParams.alpha = zeros(nNeurons, 1) + 2;   % (fixed) Aspect ratio - controls sharpness. Keep this fixed. Reducing the value makes the changes very rapid towards the end which we probably don't want.
tuningParams.b     = zeros(nNeurons, 1) + 2;   % (fixed maybe/variable - (0.5, some max - 3, 4 ...)) Control this - Control sharpness + range of the neuron. Set it to 2 for these simulations
tuningParams.q     = zeros(nNeurons, 1) + 2;   % 1 or 2?? (variable) Set it to some constant. Controls the sharpness and amplitude of peak FR.
tuningParams.w     = zeros(nNeurons, 1) + 0;   % 0 or 1?? Doesn't matter since untuned component is zero. Weight of untuned filter. Set it to 0. (fixed) Doesn't matter what is val is becz untuned filter amp is zero.
tuningParams.e1    = zeros(nNeurons, 1);       % Stimulus independent spontaneous discharge (variable) Controls dynamic range.
tuningParams.e2    = zeros(nNeurons, 1);       % Stimulus dependent spontaneous discharge (variable) Controls dynamic range.
tuningParams.gam   = zeros(nNeurons, 1);       % Controls response amplitude (variable) Controls dynamic range.
tuningParams.beta  = zeros(nNeurons, 1) + 0;   % Stimulus independent constant (variable) Controls dynamic range.
tuningParams.UNTUNED_FILTER_AMPL = 0;          % (fixed) Untuned filter not needed.
tuningParams.eps   = 0;                        % (this probably needs to be sampled every trial??) Not right place to update here Normalization noise sampled from some distribution with sigma_g standard deviation
tuningParams.normalizationNoiseParam = 0;

values = 1:8;
% probabilities (peak at 2)
% p = [0.1 0.3 0.2 0.15 0.1 0.07 0.05 0.03];
p = [0.1 0.3 0.2 0.15 0.1 0.07 0.05 0.03];
p = p / sum(p); % normalize
samples = randsample(values,100,true,p);

% Overriding variable tuning parameters
% These parameters remain same for all the trials i.e. it doesn not change
% form trial by trial
tuningParams.alpha(:)   = 2; %2; % Varying this creates non-uniform tuning curves 1 + 5*rand(nNeurons, 1);            % 1 + 5*rand(1, nNeurons);  Aspect ratio uniformly sampled from 0 - 5
tuningParams.b(:)       = 2; %randsample(values,nNeurons,true,p); % Derivative order - Just set it 2 (Zoey's paper). Non integer vaalues can give imaginary values. More like log uniform 0.0125 + 8*rand(1, nNeurons); 
tuningParams.q(:)       = exp( log(1.8)*rand(nNeurons, 1) );  % Transduction (this is not uniformly distributed) exp( log(1.8)*rand(1, nNeurons) )
tuningParams.beta(:)    = lognrnd(2.5, 0.5, nNeurons, 1);     % lognrnd(2.5, 0.5, 1, nNeurons);     % Normalization constant lognrnd(mu, sigma, 1, nneurons)

tuningParams.e1(:)        = lognrnd(0.8, 0.6, nNeurons, 1);       % Range: 0 - 10 ips, Might not be correct initilization
% tuningParams.e2(:)      = lognrnd(0.8, 0.6, 1, nNeurons);
tuningParams.gam(:)       = 4000; % This is good parameter values
neuronsPrefOrientation(:) = pi * rand(nNeurons, 1);        % Randomly choose neurons preferred orientation from 0 to pi (not directional selective)

% 1, 3
% tuningParams.alpha(:)   = 1 + 1*rand(1, nNeurons);            % Aspect ratio uniformly sampled from 0 - 5
% tuningParams.b(:)       = 2;                                  % Derivative order - Just set it 2 (Zoey's paper). Non integer vaalues can give imaginary values. More like log uniform 0.0125 + 8*rand(1, nNeurons); 
% tuningParams.q(:)       = exp( log(1)*rand(1, nNeurons) );    % Transduction (this is not uniformly distributed)
% tuningParams.beta(:)    = lognrnd(2.5, 0.5, 1, nNeurons); %lognrnd(2.5, 0.5, 1, nNeurons);     % Normalization constant lognrnd(mu, sigma, 1, nneurons)
% 
% tuningParams.e1(:)   = lognrnd(1, 0.8, 1, nNeurons);       % Range: 0 - 10 ips, Might not be correct initilization
% % tuningParams.e2(:) = lognrnd(1, 0.8, 1, nNeurons);
% tuningParams.gam(:)  = 500;
% neuronsPrefOrientation(:) = pi * rand(1, nNeurons);        % Randomly choose neurons preferred orientation from 0 to pi (not directional selective)


% ----------------------------------
% Computing stimulus response begins
% ----------------------------------

% STEP1: 
% Does this for reach trial in a loop
% Loop over every trial instead of ending one single vector. Normalization
% noise needs to be sampled every trial which is part of the tuning params.
% stimParams.contrastLevel = 1;
% stimParams.spreadLevel = deg2rad(10);
% stimParams.stimOri = deg2rad(90); % does the peak stim orientation needs to be noisy i.e. center shifts from trial to trial (probably yes)?
% 
% [tuningFnData] = getOriTunedStimRespFuntion( ...
%     neuronsPrefOrientation, tuningParams, stimParams);

stimOris    = linspace(0, pi, 201);
contrasts   = [0.01 0.02]; %   0.05 0.01 0.05 0.1 It saturated after 0.1
dispersion  = [3 30];  % 30
tuningFn    = zeros(nNeurons, numel(contrasts), numel(dispersion), numel(stimOris));

tuningFnDatas          = cell(1, 2);
stimFltResp            = zeros(numel(contrasts), numel(dispersion), numel(stimOris), nNeurons);
unnormalizedTuningFn   = zeros(numel(contrasts), numel(dispersion), numel(stimOris), nNeurons, 501);
stimProfile            = zeros(numel(contrasts), numel(dispersion), 501);

for i=1:numel(stimOris)
    for c=1:numel(contrasts)
        for d=1:numel(dispersion)
            
            stimParams.contrastLevel = contrasts(c); %1e-11; % 1e-11 Keep this small
            stimParams.spreadLevel = deg2rad(dispersion(d));
            stimParams.stimOri = stimOris(i);
            
            [tuningFnData] = getOriTunedStimRespFunction( ...
                neuronsPrefOrientation, tuningParams, stimParams);
            
            tuningFn(:, c, d, i) = tuningFnData.FR;
            
            tuningFnDatas{d} = tuningFnData;
            stimFltResp(c, d, i, :) = tuningFnData.stimfltResp;
            unnormalizedTuningFn(c, d, i, :, :) = tuningFnData.unnormalizedResp;
            stimProfile(c, d, :) = tuningFnData.stimProfile;
        end
    end   
end

%% Plots

for c=1:numel(contrasts)
    for d=1:numel(dispersion)
        tf        = squeeze( tuningFn(:, c, d, :) );
        OSI       = abs(tf*(exp(1i*(2*stimOris))'))./sum(abs(tf), 2);
        peakSpkRt = squeeze(max(tf, [], 2));
        
        figure

        subplot(2, 2, 1)
        hold on
        histogram(OSI, 10, 'Normalization', 'probability') % 'Normalization', 'probability'
        xlabel('Orientation Selectivity (OSI)')
        ylabel('Proportion')
        xlim([0, 1])
        ylim([0, 0.5])
        hold off
        title(sprintf("C: %.2f, D: %d", contrasts(c), dispersion(d)))
        
        % Plot all tuning curves
        subplot(2, 2, 2)
        hold on
        title("Tuning curves" + newline + " (all neurons)")
        for k = 1:nNeurons
            plot(rad2deg(stimOris), squeeze(tuningFn(k, c, d, :)))
        end
        xlabel('Orientation (deg)')
        ylabel('ips')
        % ylim([0 1])
        % xlim([0, 180])
        hold off
        
        subplot(2, 2, 3)
        hold on
        histogram(peakSpkRt, 50, 'Normalization', 'probability')
        xlabel('Peak IPS (evoked)')
        ylabel('Probability')
        % xlim([0 40])
        hold off
        
        subplot(2, 2, 4)
        hold on
        histogram(tuningParams.e1, 20, 'Normalization', 'probability')
        xlabel('Baseline FR (spontaneous)')
        ylabel('Probability')
        hold off

    end
end   

%% Stim filter response
figure

for c=1:numel(contrasts)
    for d=1:numel(dispersion)
        sResp    = squeeze( stimFltResp(c, d, :, :) );
        sumsResp = sum(sResp, 2);
        
        subplot(2, 2, 1)
        hold on
        histogram(sumsResp(:), DisplayName=sprintf("C: %.2f, D: %d", contrasts(c), dispersion(d)))
        xlabel("Filter stim response (aggregate over all the neurons)")
        ylabel("count")
        title("Resp at each orientation")
        legend
        set(gca, 'XScale', 'log')
        hold off
        
        subplot(2, 2, 2)
        hold on
        histogram(sResp(:), DisplayName=sprintf("C: %.2f, D: %d", contrasts(c), dispersion(d)))
        xlabel("Filter stim response")
        ylabel("count")
        title("Resp at each orienation for each neuron")
        legend
        set(gca, 'XScale', 'log')
        hold off
    end
end

for c=1:numel(contrasts)
    idx = 1;
    sResp    = squeeze( stimFltResp(c, idx, :, :) );
    sumsResp = sum(sResp, 2);
    
    subplot(2, 2, 3)
    hold on
    histogram(sumsResp(:), DisplayName=sprintf("C: %.2f", contrasts(c)))
    xlabel("Filter stim response (aggregate over all the neurons)")
    ylabel("count")
    title(sprintf("Aggregate resp at dispersion %d", dispersion(idx)))
    legend
    set(gca, 'XScale', 'log')
    hold off
end

for d=1:numel(dispersion)
    idx = 1;
    sResp    = squeeze( stimFltResp(idx, d, :, :) );
    sumsResp = sum(sResp, 2);
    
    subplot(2, 2, 4)
    hold on
    histogram(sumsResp(:), DisplayName=sprintf("D: %.2f", dispersion(d)))
    xlabel("Filter stim response (aggregate over all the neurons)")
    ylabel("count")
    title(sprintf("Aggregate resp at contrast: %.2f ", contrasts(idx)))
    legend
    set(gca, 'XScale', 'log')
    hold off
end


%% Parameters
figure

subplot(3, 3, 1)
histogram(tuningParams.d(:))
ylim([0 300])
xlim([ min(tuningParams.d(:)) - 1, max(tuningParams.d(:)) + 1 ])
title("Direction selectivity (d)")

subplot(3, 3, 2)
histogram(tuningParams.alpha(:))
ylim([0 300])
title("Aspect ratio (alpha)")

subplot(3, 3, 3)
histogram(tuningParams.b(:))
ylim([0 300])
xlim([ min(tuningParams.b(:)) - 1, max(tuningParams.b(:)) + 1 ])
title("Derivative order (b)")

subplot(3, 3, 4)
histogram(tuningParams.q(:))
ylim([0 300])
title("Transduction (q)")

subplot(3, 3, 5)
histogram(tuningParams.w(:))
ylim([0 300])
xlim([ min(tuningParams.w(:)) - 1, max(tuningParams.w(:)) + 1 ])
title("Untuned component (w)")

subplot(3, 3, 6)
histogram(tuningParams.e1(:))
ylim([0 300])
title("Spontanoues discharge (e1)")

subplot(3, 3, 7)
histogram(tuningParams.e2(:))
ylim([0 300])
xlim([ min(tuningParams.e2(:)) - 1, max(tuningParams.e2(:)) + 1 ])
title("Spontanoues discharge (e2)")

subplot(3, 3, 8)
histogram(tuningParams.gam(:))
ylim([0 300])
title("Response amplitude (gamma)")

subplot(3, 3, 8)
histogram(tuningParams.beta(:))
ylim([0 300])
title("Stim independent cont (Beta)")


%
figure
% Plot all tuning curves
hold on
title("Tuning curves" + newline + " (all neurons)")
for i = 1:20 %nNeurons
    plot(rad2deg(tuningFnData.oriSpace), tuningFnData.unnormalizedResp(i, :))
end
xlabel('Orientation (deg)')
ylabel('ips')
% xlim([0, 180])
hold off


% TODO: mean vs variance realtion


%% Tuning params with contrast and dispersion manipulation

%% Contrast tuning functions

nrnIdxes = randi([1 100], 1, 14);

%%

for d=1:numel(dispersion)
     figure
     for n=1:numel(nrnIdxes)
        subplot(4, 4, n)
    
        for c=1:numel(contrasts)
            % d = 1;
            tf = squeeze( tuningFn(n, c, d, :) );
            hold on
            plot(rad2deg(stimOris), tf, DisplayName=sprintf("C: %.2f", contrasts(c)), LineWidth=1.5)
            hold off
            xlabel("oreintation")
            ylabel("FR")
            title(sprintf("Dispersion: %d", dispersion(d)))
            legend
        end
    end

end

%% Dispersion tuning functions

for c = 1:numel(contrasts)
    figure
    for n=1:numel(nrnIdxes)
        subplot(4, 4, n)
    
        for d=1:numel(dispersion)
            % c = 1;
            tf = squeeze( tuningFn(n, c, d, :) );
            hold on
            plot(rad2deg(stimOris), tf, DisplayName=sprintf("D: %d", dispersion(d)), LineWidth=1.5)
            hold off
            xlabel("oreintation")
            ylabel("FR")
            title(sprintf("Contrast: %.2f", contrasts(c)))
            legend
        end
    end
end

%% Tuning function and stim response

c = 1;
ori = tuningFnData.oriSpace;
x1  = squeeze( unnormalizedTuningFn(c, 1, 1, 1, :) );
x2  = squeeze( unnormalizedTuningFn(c, 2, 1, 1, :) );
p1  = squeeze( stimProfile(c, 1, :) );
p2  = squeeze( stimProfile(c, 2, :) );

figure
hold on
plot(rad2deg(ori), x1, DisplayName="Tuning fn")
plot(rad2deg(ori), p1, DisplayName="Stim profile")
hold off
xlabel("orientation")
ylim([0 3])
legend
title(sprintf("dispersion %d", dispersion(1)))

figure
hold on
plot(rad2deg(ori), x2, DisplayName="Tuning fn")
plot(rad2deg(ori), p2, DisplayName="Stim profile")
hold off
xlabel("orientation")
ylim([0 3])
legend
title(sprintf("dispersion %d", dispersion(2)))

