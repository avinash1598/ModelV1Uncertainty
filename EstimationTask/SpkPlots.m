clear all
% close all
clc

rng('shuffle');

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/ModelV1Uncertainty/Model/Scripts/')

% TIP: Keep all possible 1D variables to row vector
% Note: Time step should be set carefully. It should be pretty small
% relative to the spike rate of single neuron so that the probability of
% firing does not shoot up.
% TODO: add a check to make sure firing rate is not so big compared to the
% time window.

% ----------------------------------
% Params
% ----------------------------------
nNeurons     = 200;   % Count of neurons
stimDuration = 20;     % Stimulus duration set to 1 seconds
timeStep     = 0.01; % 0.001s (1ms) - Step size of time bins used for binning stimulus duration 
% ntrials    = 100;   % Offloaded to generateStim function. No of stimulus to simulate the model with

% Stimulus parameters (angles in radians)
contrasts               = [0.01 0.05]; % 0.01
spreads                 = [3 30];
uniqStimOris            = 0:30:180;
uniqStimOris            = deg2rad(uniqStimOris');
stimParam.numStim       = numel(uniqStimOris);                                % Number of unique stimuli
stimParam.countPerStim  = 100; 
ntrials                 = stimParam.numStim * stimParam.countPerStim * numel(contrasts) * numel(spreads);  % Total number of trials

[c, s, oris] = ndgrid(contrasts, spreads, uniqStimOris);
combinations = [c(:), s(:) oris(:)];
trlMatrix    = repmat(combinations, [stimParam.countPerStim 1]);

%%
% Add random noise to the stimulus orientation
% Stim noise can itself cause gain variability
trlContrastVector = trlMatrix(:, 1);
trlSpreadVector   = trlMatrix(:, 2);
trlStimVector     = trlMatrix(:, 3);  % Vector of stimuli

stimNoise       = 0 + 0 * randn(ntrials, 1);                    % What is std dev here? 5.73 degrees 0.1 * randn(ntrials, 1); I guess even noise in the stimulus can cause modulated poisson process  
noisyStimVector = trlStimVector + stimNoise;                         % Noisy stimulus vector


% Shuffle stim vector - this is important for studying effect of top down
% effect
shuffleIdx        = randperm(ntrials);
trlContrastVector = trlContrastVector(shuffleIdx);
trlSpreadVector   = trlSpreadVector(shuffleIdx);
trlStimVector     = trlStimVector(shuffleIdx);
noisyStimVector   = noisyStimVector(shuffleIdx); 

neuronsPrefOrientation = zeros(nNeurons, 1);
timeBins               = 0:timeStep:stimDuration; 
stimRespProfile        = 1 + zeros(1, numel(timeBins));
gainVector             = 1 + zeros(nNeurons, ntrials); % constant gain - NO gain modulation

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
tuningParams.normalizationNoiseParam = 5; %230;      % Std dev in radians sigma_n

values = 1:8;
% probabilities (peak at 2)
% p = [0.1 0.3 0.2 0.15 0.1 0.07 0.05 0.03];
p = [0.1 0.3 0.2 0.15 0.1 0.07 0.05 0.03];
p = p / sum(p); % normalize
samples = randsample(values,100,true,p);

% Overriding variable tuning parameters
% These parameters remain same for all the trials i.e. it doesn not change
% form trial by trial
tuningParams.alpha(:)   = 2; %1 + 5*rand(nNeurons, 1);            % 1 + 5*rand(1, nNeurons);  Aspect ratio uniformly sampled from 0 - 5
tuningParams.b(:)       = 2; %randsample(values,nNeurons,true,p); % Derivative order - Just set it 2 (Zoey's paper). Non integer vaalues can give imaginary values. More like log uniform 0.0125 + 8*rand(1, nNeurons); 
tuningParams.q(:)       = exp( log(1.8)*rand(nNeurons, 1) );  % Transduction (this is not uniformly distributed) exp( log(1.8)*rand(1, nNeurons) )
tuningParams.beta(:)    = lognrnd(2.5, 0.5, nNeurons, 1);     % lognrnd(2.5, 0.5, 1, nNeurons);     % Normalization constant lognrnd(mu, sigma, 1, nneurons)

tuningParams.e1(:)        = lognrnd(0.8, 0.6, nNeurons, 1);       % Range: 0 - 10 ips, Might not be correct initilization
% tuningParams.e2(:)      = lognrnd(0.8, 0.6, 1, nNeurons);
tuningParams.gam(:)       = 4000;
neuronsPrefOrientation(:) = pi * rand(nNeurons, 1);        % Randomly choose neurons preferred orientation from 0 to pi (not directional selective)


%% Generate spike data
% Structures to store final neuron spikes
% Preallocate result and spike response matrices
% firingRates          = squeeze( tuningFn(:, 1, 1, :) );
trialDecisions       = zeros(ntrials, 1);
neuronSpikeResponses = false(ntrials, nNeurons, length(timeBins)); % Creating a logical matrix to save memory

for trialIDx = 1:ntrials
    if mod(trialIDx, stimParam.countPerStim) == 0
        disp(trialIDx)
    end
    
    % Get tuning params every trial
    stimParams.contrastLevel = trlContrastVector(trialIDx); %0.01 0.04 0.1 Use these two values. Saturation happens pretty quickly
    stimParams.spreadLevel   = deg2rad(trlSpreadVector(trialIDx));
    stimParams.stimOri       = noisyStimVector(trialIDx);
    
    [tuningFnData] = getOriTunedStimRespFuntion( ...
            neuronsPrefOrientation, tuningParams, stimParams);
    firingRates = tuningFnData.FR; % Firing rate for this trial
    
    % STEP 2: Compute stimulus response for each neuron over time
    % This multiplies the firing rates with a time-dependent stimulus response profile
    % Output: 
    %  - trlStimResponse: response of each neuron over time for each trial (nNeurons x nTimeBins)
    % trlStimResponse = firingRates(:, trialIDx).*stimRespProfile;
    trlStimResponse = firingRates.*stimRespProfile;
    
    % Modulate gain of current trial based on previous trial
    % nNeurons x No time bins
    trlGainVector = squeeze(repmat(gainVector(:, trialIDx)', [1, 1, length(timeBins)])); % Extract gain vector for this trial for each timebin
    
    % STEP 3: Generate modulated Poisson spikes for each trial
    % Output:
    %  - spikes: spike trains for each neuron in the trial
    %  - modStimResponse: modified stimulus response after gain modulation
    params = struct();
    params.timeStep = timeStep;
    params.timeBins = timeBins;
    params.nNeurons = nNeurons;
    [spikes, modStimResponse] = generateModulatedPoissonSpikes(trlStimResponse, ...
        trlGainVector, params);
    
%     % STEP 4: Decode the stimulus orientation based on the spike trains
%     % Output:
%     %  - thetaMLE: maximum likelihood estimate of stimulus orientation based on spikes
%     %  - decodingError: error between decoded orientation and actual stimulus
%     params.stimDuration = stimDuration;
%     thetaMLE = decodeOrientationFromSpikes(spikes, ...
%         neuronsPrefOrientation, params, tuningParams);
%     decodingError = thetaMLE - noisyStimVector(trialIDx);
%     
%     % Decision: CW (-1) or CCW (1) based on decoded orientation
%     decision = (thetaMLE(end) > pi/2)*(-1) + (thetaMLE(end) <= pi/2)*(1);
    
    % Store spike trains and decision results
    neuronSpikeResponses(trialIDx, :, :) = logical(spikes);  % Store spikes
    % trialDecisions(trialIDx) = decision;  % Store decision result (CW or CCW)
end


%% Figures

figure
subplot(2, 3, 1)
hold on
title("Single trial response" + newline + "all neurons")
for nIDx = 1:nNeurons
    plot(squeeze(timeBins), squeeze(modStimResponse(nIDx, :)));
end
axis square
xlabel("Time (s)")
ylabel("IPS")
hold off

subplot(2, 3, 2)
hold on
title("Single trial spikes" + newline + "all neurons")
imagesc(squeeze(spikes(:, :))), colormap(flipud('gray'))
axis square
box off, axis off
xlabel("Time (s)")
ylabel("Spikes")
hold off

%%
% Plot aggregate
% Plot contrasts
% Plot spreads

% Sanity check: Mean vs variance plot
% TODO: do it grouped by single orientation??
numIntervals = 1;
intervalSize = floor(length(timeBins) / numIntervals);  % 20 columns per interval
meanSpkCnt = zeros(nNeurons, numel(contrasts), numel(spreads), numel(uniqStimOris), numIntervals);
varSpkCnt  = zeros(nNeurons, numel(contrasts), numel(spreads), numel(uniqStimOris), numIntervals);

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        for oIdx = 1:numel(uniqStimOris)
            stimVal  = uniqStimOris(oIdx);
            trlIdxes = ( trlStimVector ==  stimVal ) & (trlContrastVector == contrasts(cIdx)) & (trlSpreadVector == spreads(sIdx));
            
            for nIdx = 1:nNeurons
                spkForThisNrn = squeeze( neuronSpikeResponses(trlIdxes, nIdx, :) );
                
                for i = 1:numIntervals
                    startCol = (i-1) * intervalSize + 1;
                    endCol = min(i * intervalSize, length(timeBins));  % Handle the last interval
                    
                    % mean spk rate
                    spkCntInThisInterval = sum( spkForThisNrn(:, startCol:endCol), 2 );
                    intervalData = mean( spkCntInThisInterval ); 
                    
                    meanSpkCnt(nIdx, cIdx, sIdx, oIdx, i) = mean( spkCntInThisInterval );
                    varSpkCnt(nIdx, cIdx, sIdx, oIdx, i)  = var( spkCntInThisInterval );
                end
            end
        end
    end
end

figure

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        hold on
        m = meanSpkCnt(:, cIdx, sIdx, :); m = m(:);
        v = varSpkCnt(:, cIdx, sIdx, :); v = v(:);
        
        scatter(m, v, DisplayName=sprintf("C: %.2f, D: %d", contrasts(cIdx), spreads(sIdx)))
        hold on
        plot([min(m) max(m)], [min(m) max(m)], 'k--', HandleVisibility='off') 
        hold off
        xlabel("Mean spk cnt")
        ylabel("Var spk cnt")
        title("Conditioned on stimulus orientation")
        legend
        set(gca, 'XScale', 'log', 'YScale', 'log')
        hold off
    end
end

figure
for cIdx = 1:numel(contrasts)
    subplot(2, 2, cIdx)
    hold on

    for sIdx = 1:numel(spreads)
        
        m = meanSpkCnt(:, cIdx, sIdx, :); m = m(:);
        v = varSpkCnt(:, cIdx, sIdx, :); v = v(:);
        
        scatter(m, v, DisplayName=sprintf("D: %d", spreads(sIdx)))
        plot([min(m) max(m)], [min(m) max(m)], 'k--', HandleVisibility='off') 
        
    end

    xlabel("Mean spk cnt")
    ylabel("Var spk cnt")
    title(sprintf("Contrast: %.2f", contrasts(cIdx)))
    legend
    set(gca, 'XScale', 'log', 'YScale', 'log')
    hold off
end

for sIdx = 1:numel(spreads)
    subplot(2, 2, numel(contrasts) + sIdx)
    hold on
    
    for cIdx = 1:numel(contrasts)
        
        m = meanSpkCnt(:, cIdx, sIdx, :); m = m(:);
        v = varSpkCnt(:, cIdx, sIdx, :); v = v(:);
        
        scatter(m, v, DisplayName=sprintf("C: %.2f", contrasts(cIdx)))
        plot([min(m) max(m)], [min(m) max(m)], 'k--', HandleVisibility='off') 
        
    end
    
    xlabel("Mean spk cnt")
    ylabel("Var spk cnt")
    title(sprintf("Spread: %d", spreads(sIdx)))
    legend
    set(gca, 'XScale', 'log', 'YScale', 'log')
    hold off
end

% numIntervals = 30;
% intervalSize = floor(length(timeBins) / numIntervals);  % 20 columns per interval
% meanSpkCnt = zeros(nNeurons, numel(uniqStimOris), numIntervals);
% varSpkCnt  = zeros(nNeurons, numel(uniqStimOris), numIntervals);
% 
% for nIdx = 1:nNeurons
%     spkForThisNrn = squeeze( neuronSpikeResponses(:, nIdx, :) );
%     
%     for i = 1:numIntervals
%         startCol = (i-1) * intervalSize + 1;
%         endCol = min(i * intervalSize, length(timeBins));  % Handle the last interval
%         
%         % mean spk rate
%         spkCntInThisInterval = sum( spkForThisNrn(:, startCol:endCol), 2 );
%         intervalData = mean( spkCntInThisInterval ); 
%         
%         meanSpkCnt(nIdx, i) = mean( spkCntInThisInterval );
%         varSpkCnt(nIdx, i)  = var( spkCntInThisInterval );
%     end
% end
% 
% 
% subplot(2, 2, 2)
% scatter(meanSpkCnt(:), varSpkCnt(:))
% hold on
% plot([0 max(meanSpkCnt(:))], [0 max(meanSpkCnt(:))], 'k--') 
% hold off
% xlabel("Mean spk cnt")
% ylabel("Var spk cnt")
% title("Each point - spikes from multiple orientation values")


% Normalization noise -
% I think normalization noise should be constant and should not depend upon
% uncertainty level.
% 1. normalization noise - is it same for all the neurons
% 1. probably needs to be uncertainty level dependent
% 2. It is picked from gaussian distribution - need to make sure the net
% summation in the denominator is not negative. finding this value might be
% hand wavey since normalization pool summation depends upon the stimulus
% itself
% 3. Other option to introduce normalization noise is to pick the stimulus
% orientation from a distribution dependent upon uncertainty level
% 4. not sure if this is right from  decoding perspective