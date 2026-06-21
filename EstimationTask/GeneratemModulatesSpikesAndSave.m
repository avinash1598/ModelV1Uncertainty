clear all
close all
clc

rng('shuffle');

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/ModelV1Uncertainty/Model/Scripts/')

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
nNeurons     = 200;   % Count of neurons
stimDuration = 2;     % Stimulus duration set to 1 seconds
timeStep     = 0.001; % 0.001s (1ms) - Step size of time bins used for binning stimulus duration 
% ntrials    = 100;   % Offloaded to generateStim function. No of stimulus to simulate the model with

% Stimulus parameters (angles in radians)
contrasts               = [0.01 0.05]; % 0.01
spreads                 = [3 30];
% sigma_c                 = [10 0.5]; % sigma contribution from contrast
% sigma_s                 = [5 15]; % sigma contribution from spread
noise                   = deg2rad( [10 3 51 16] ); % 30 5 31 21 For now just hardcoding - (0.01, 3) (0.05, 3) (0.01, 30) (0.05, 30) 
uniqStimOris            = 0:10:180;
uniqStimOris            = deg2rad(uniqStimOris');
stimParam.numStim       = numel(uniqStimOris);                                % Number of unique stimuli
stimParam.countPerStim  = 100; 
ntrials                 = stimParam.numStim * stimParam.countPerStim * numel(contrasts) * numel(spreads);  % Total number of trials

noise        = repmat(noise', [stimParam.numStim 1]);
[c, s, oris] = ndgrid(contrasts, spreads, uniqStimOris);
combinations = [c(:), s(:) oris(:) noise(:)];
trlMatrix    = repmat(combinations, [stimParam.countPerStim 1]);

[c, s, oris] = ndgrid(contrasts, spreads, uniqStimOris);
combinations = [c(:), s(:) oris(:)];

% % Set noise level for each conditions
% noise = [];
% for idx = 1:size(combinations, 1)
%     sig1 = sigma_c(contrasts == c(idx));
%     sig2 = sigma_s(spreads == s(idx));
%     
%     noise = [noise sqrt(sig1^2 + sig2^2)];   
% end
% 
% noise        = deg2rad(noise);
% combinations = [combinations noise(:)];
% trlMatrix    = repmat(combinations, [stimParam.countPerStim 1]);

%% TODO: do this in a loop to make sure its right

% Add random noise to the stimulus orientation
% Stim noise can itself cause gain variability
% Don't forget to shuffle
trlContrastVector = trlMatrix(:, 1);
trlSpreadVector   = trlMatrix(:, 2);
trlStimVector     = trlMatrix(:, 3);  % Vector of stimuli
stimNoiseVector   = 0 + trlMatrix(:, 4) .* randn(ntrials, 1);   % This is just changin mean value. change the width of stim response profile as well.
noisyStimVector   = trlStimVector + stimNoiseVector;                         % Noisy stimulus vector

% Shuffle stim vector - this is important for studying effect of top down
% effect
shuffleIdx        = randperm(ntrials);
trlContrastVector = trlContrastVector(shuffleIdx);
trlSpreadVector   = trlSpreadVector(shuffleIdx);
trlStimVector     = trlStimVector(shuffleIdx);
stimNoiseVector   = stimNoiseVector(shuffleIdx);
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
% tuningParams.normalizationNoiseParam = 0; %230;      % Std dev in radians sigma_n

% values = 1:8;
% % probabilities (peak at 2)
% % p = [0.1 0.3 0.2 0.15 0.1 0.07 0.05 0.03];
% p = [0.1 0.3 0.2 0.15 0.1 0.07 0.05 0.03];
% p = p / sum(p); % normalize
% samples = randsample(values,100,true,p);

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

%% Tuning function
data = {};

data.stimDuration = stimDuration;    
data.timeStep     = timeStep; 
data.nNeurons     = nNeurons;
data.contrasts    = contrasts;
data.spreads      = spreads;
data.uniqStimOris = uniqStimOris;
data.tuningParams = tuningParams;
data.neuronsPrefOrientation = neuronsPrefOrientation;

% No need to save tuning functions. Just save the parameters.

% for cIdx = 1:numel(contrasts)
%     for sIdx = 1:numel(spreads)
% 
%         tuningFns = [];
%         
%         for oIdx = 1:numel(uniqStimOris)
%             
%             % Get tuning params every trial
%             stimParams.contrastLevel = contrasts(cIdx);
%             stimParams.spreadLevel   = deg2rad(spreads(sIdx));
%             stimParams.stimOri       = uniqStimOris(oIdx);
%             
%             tFn = getOriTunedStimRespFunction( ...
%                     neuronsPrefOrientation, tuningParams, stimParams);
%             firingRates = tFn.FR; % Firing rate for this trial
%             
%             tuningFns = [tuningFns firingRates];
%         end
% 
%         key1 = sprintf('c_%g', contrasts(cIdx));
%         key2 = sprintf('s_%g', spreads(sIdx));
%         key1 = matlab.lang.makeValidName(key1);
%         key2 = matlab.lang.makeValidName(key2);
%         data.tuningFns.(key1).(key2) = tuningFns;
%     end
% end

% save('tuningData.mat', 'data');

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
    
    tFn = getOriTunedStimRespFunction( ...
            neuronsPrefOrientation, tuningParams, stimParams);
    firingRates = tFn.FR; % Firing rate for this trial
    
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


%%
% calculate gain for each neuron and for each condition
% load("SpkData.mat")

% meanSpkCnt = dataToSave.meanSpkCnt;
% varSpkCnt  = dataToSave.varSpkCnt;


for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)

        gainVals = [];

        for n=1:nNeurons
            muSpkCnt     = meanSpkCnt(n, cIdx, sIdx, :);
            sigma2SpkCnt = varSpkCnt(n, cIdx, sIdx, :);

            x = muSpkCnt(:); 
            y = sigma2SpkCnt(:);

            % Define custom model (edit as needed)
            ft = fittype('x + sigma_g^2*(x)^2', ...
                         'independent','x','coefficients',{'sigma_g'});
            
            opts = fitoptions(ft);
            opts.StartPoint = [rand];     % initial guess
            opts.Lower = [0];          % e.g., sigma_g >= 0
            % opts.Upper = [10];       % optional upper bound
            
            [curve, ~] = fit(x, y, ft, opts);
            coeffs = coeffvalues(curve);
            
            gainVals = [gainVals coeffs];

        end

        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);
        data.nrnGains.(key1).(key2) = gainVals;
    end
end

data.meanSpkCnt = meanSpkCnt;
data.varSpkCnt = varSpkCnt;

save('tuningFnData.mat', 'data');


