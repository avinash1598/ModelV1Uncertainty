clear all
close all
clc

rng('shuffle');

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/ModelV1Uncertainty/Model/Scripts/')

tuningFnData = load('tuningFnData.mat');

% IMP: For gain variability to increase with increasing dispersion, sum of
% normalization signal should decrease. But that does not necessarily
% happens with dispersion. It is always true for contrast though.

% TIP: Keep all possible 1D variables to row vector
% Note: Time step should be set carefully. It should be pretty small
% relative to the spike rate of single neuron so that the probability of
% firing does not shoot up.
% TODO: add a check to make sure firing rate is not so big compared to the
% time window.

% ----------------------------------
% Params
% ----------------------------------
nNeurons     = tuningFnData.data.nNeurons;      % Count of neurons
stimDuration = tuningFnData.data.stimDuration;  % Stimulus duration set to 1 seconds
timeStep     = tuningFnData.data.timeStep;      % 0.001s (1ms) - Step size of time bins used for binning stimulus duration 

% Stimulus parameters (angles in radians)
contrasts               = tuningFnData.data.contrasts; %[0.01 0.05]; % 0.01
spreads                 = tuningFnData.data.spreads;   %[3 30];
% noise                   = deg2rad( [10 3 51 16] ); % 30 5 31 21 For now just hardcoding - (0.01, 3) (0.05, 3) (0.01, 30) (0.05, 30) 
varGains                = [];
uniqStimOris            = tuningFnData.data.uniqStimOris; % 0:10:180;
%uniqStimOris            = deg2rad(uniqStimOris');
stimParam.numStim       = numel(uniqStimOris);                                % Number of unique stimuli
stimParam.countPerStim  = 10; % 100 
ntrials                 = stimParam.numStim * stimParam.countPerStim * numel(contrasts) * numel(spreads);  % Total number of trials

% noise        = repmat(noise', [stimParam.numStim 1]);
[c, s, oris] = ndgrid(contrasts, spreads, uniqStimOris);
combinations = [c(:), s(:) oris(:) noise(:)];
trlMatrix    = repmat(combinations, [stimParam.countPerStim 1]);

%% TODO: do this in a loop to make sure its right

% Add random noise to the stimulus orientation
% Stim noise can itself cause gain variability
% Don't forget to shuffle
trlContrastVector = trlMatrix(:, 1);
trlSpreadVector   = trlMatrix(:, 2);
trlStimVector     = trlMatrix(:, 3);  % Vector of stimuli
stimNoiseVector   = 0 + trlMatrix(:, 4) .* randn(ntrials, 1);      % This is just changin mean value. change the width of stim response profile as well.
noisyStimVector   = trlStimVector + stimNoiseVector;               % Noisy stimulus vector
noisyStimVector   = deg2rad( mod(rad2deg(noisyStimVector), 180) ); % Wrap between 0 and 180

% Shuffle stim vector - this is important for studying effect of top down
% effect
shuffleIdx        = randperm(ntrials);
trlContrastVector = trlContrastVector(shuffleIdx);
trlSpreadVector   = trlSpreadVector(shuffleIdx);
trlStimVector     = trlStimVector(shuffleIdx);
stimNoiseVector   = stimNoiseVector(shuffleIdx);
noisyStimVector   = noisyStimVector(shuffleIdx); 

neuronsPrefOrientation = tuningFnData.data.neuronsPrefOrientation; % zeros(nNeurons, 1);
timeBins               = 0:timeStep:stimDuration; 
stimRespProfile        = 1 + zeros(1, numel(timeBins));
gainVector             = 1 + zeros(nNeurons, ntrials); % constant gain - NO gain modulation
tuningParams           = tuningFnData.data.tuningParams;

%% Update tuning function data
tuningFnOriSpace = linspace(0, pi, 361);

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)

        tuningFns = [];
        
        for oIdx = 1:numel(tuningFnOriSpace)
            
            % Get tuning params every trial
            stimParams.contrastLevel = contrasts(cIdx);
            stimParams.spreadLevel   = deg2rad(spreads(sIdx));
            stimParams.stimOri       = tuningFnOriSpace(oIdx);
            
            tFn = getOriTunedStimRespFunction( ...
                    neuronsPrefOrientation, tuningParams, stimParams);
            firingRates = tFn.FR; % Firing rate for this trial
            
            tuningFns = [tuningFns firingRates];
        end

        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);
        tuningFnData.data.tuningFns.(key1).(key2) = tuningFns;
    end
end

% Update tuning Fns
tuningFns = tuningFnData.data.tuningFns;
nrnGains  = tuningFnData.data.nrnGains;

%% Generate spike data
% Structures to store final neuron spikes
% Preallocate result and spike response matrices
% firingRates          = squeeze( tuningFn(:, 1, 1, :) );
decodedThetasPossDec        = zeros(ntrials, 1);
decodedContrastsPossDec     = zeros(ntrials, 1);
decodedSpreadsPossDec       = zeros(ntrials, 1);

decodedThetasMPossDec        = zeros(ntrials, 1);
decodedContrastsMPossDec     = zeros(ntrials, 1);
decodedSpreadsMPossDec       = zeros(ntrials, 1);

for trialIDx = 1:ntrials
    if mod(trialIDx, stimParam.countPerStim) == 0
        disp(trialIDx)
    end
    
    % Get tuning params every trial
    stimParams.contrastLevel = trlContrastVector(trialIDx); %0.01 0.04 0.1 Use these two values. Saturation happens pretty quickly
    stimParams.spreadLevel   = deg2rad(trlSpreadVector(trialIDx));
    stimParams.stimOri       = noisyStimVector(trialIDx);
    
    tFn = getOriTunedStimRespFunction( ...
            neuronsPrefOrientation, tuningParams, stimParams);
    firingRates = tFn.FR; % Firing rate for this trial
    trlStimResponse = firingRates.*stimRespProfile;
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
    
    % STEP 4: Decode the stimulus orientation based on the spike trains
    % Output:
    %  - thetaMLE: maximum likelihood estimate of stimulus orientation based on spikes
    %  - decodingError: error between decoded orientation and actual stimulus
    params.stimDuration = stimDuration;
    params.contrasts    = contrasts;
    params.spreads      = spreads;
    params.uniqStimOris = tuningFnOriSpace;
    
    % Poisson decoder
    [contrastMLE, spreadMLE, thetaMLE, ~] = decodePoissonSpikes( ...
        spikes, tuningFns, params);
    decodingError = thetaMLE - trlStimVector(trialIDx); % Error wrt to actual stim ori
    
    % Decoded quantititles
    decodedThetasPossDec(trialIDx)    = thetaMLE;
    decodedContrastsPossDec(trialIDx) = contrastMLE;
    decodedSpreadsPossDec(trialIDx)   = spreadMLE;
    
    % Modulated poisson decoder
    [contrastMLE, spreadMLE, thetaMLE, ~] = decodeModulatedPoissonSpikes( ...
        spikes, tuningFns, params, nrnGains);
    
    % Decoded quantititles
    decodedThetasMPossDec(trialIDx)    = thetaMLE;
    decodedContrastsMPossDec(trialIDx) = contrastMLE;
    decodedSpreadsMPossDec(trialIDx)   = spreadMLE;
    
    % Store spike trains and decision results
    % neuronSpikeResponses(trialIDx, :, :) = logical(spikes);  % Store spikes
end

%%
figure

subplot(2, 2, 1)
scatter(rad2deg(trlStimVector), rad2deg(decodedThetasPossDec), 'filled')
xlabel("Stim theta (deg)")
ylabel("Decoded theta (deg)")

subplot(2, 2, 2)
scatter(rad2deg(noisyStimVector), rad2deg(decodedThetasPossDec), 'filled')
xlabel("Stim theta noisy (deg)")
ylabel("Decoded theta (deg)")

subplot(2, 2, 3)
scatter(trlSpreadVector + rand(ntrials, 1), decodedSpreadsPossDec + rand(ntrials, 1), 'filled')
xlabel("Stim spread (deg)")
ylabel("Decoded spread (deg)")

subplot(2, 2, 4)
scatter(trlContrastVector + 0.001*rand(ntrials, 1), decodedContrastsPossDec + 0.001*rand(ntrials, 1), 'filled')
xlabel("Stim contrast")
ylabel("Decoded contrast")


figure

subplot(2, 2, 1)
scatter(rad2deg(trlStimVector), rad2deg(decodedThetasMPossDec), 'filled')
xlabel("Stim theta (deg)")
ylabel("Decoded theta (deg)")

subplot(2, 2, 2)
scatter(rad2deg(noisyStimVector), rad2deg(decodedThetasMPossDec), 'filled')
xlabel("Stim theta noisy (deg)")
ylabel("Decoded theta (deg)")

subplot(2, 2, 3)
scatter(trlSpreadVector + rand(ntrials, 1), decodedSpreadsMPossDec + rand(ntrials, 1), 'filled')
xlabel("Stim spread (deg)")
ylabel("Decoded spread (deg)")

subplot(2, 2, 4)
scatter(trlContrastVector + 0.001*rand(ntrials, 1), decodedContrastsMPossDec + 0.001*rand(ntrials, 1), 'filled')
xlabel("Stim contrast")
ylabel("Decoded contrast")

figure
scatter(rad2deg(decodedThetasPossDec), rad2deg(decodedThetasMPossDec), 'filled')
hold on
plot([0 180], [0 180], 'k--', HandleVisibility='off') 
hold off
xlabel("Decoded theta (poisson)")
ylabel("Decoded theta (modulated poisson)")
