clear all
close all
clc

rng('shuffle');

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/ModelV1Uncertainty/EstimationTask/Scripts/')

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

% Stimulus parameters (angles in radians)
contrasts               = linspace(1e-4, 0.25, 15); %[1e-4 0.001 0.01 0.05 0.1 0.2 0.5]; % 0.01
spreads                 = linspace(1, 60, 30); %[3];
uniqStimOris            = 0:10:180;
uniqStimOris            = deg2rad(uniqStimOris');
stimParam.numStim       = numel(uniqStimOris);                                % Number of unique stimuli
stimParam.countPerStim  = 100; 
ntrials                 = stimParam.numStim * stimParam.countPerStim * numel(contrasts) * numel(spreads);  % Total number of trials

[c, s, oris] = ndgrid(contrasts, spreads, uniqStimOris);
combinations = [c(:), s(:) oris(:)];
trlMatrix    = repmat(combinations, [stimParam.countPerStim 1]);

% Add random noise to the stimulus orientation
% Stim noise can itself cause gain variability
% Don't forget to shuffle
trlContrastVector = trlMatrix(:, 1);
trlSpreadVector   = trlMatrix(:, 2);
trlStimVector     = trlMatrix(:, 3);  % Vector of stimuli
noisyStimVector   = trlStimVector;    % Noisy stimulus vector

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

%% Update tuning function data
tuningFnOriSpace = linspace(0, pi, 361);
maxFRs = zeros(numel(contrasts), numel(spreads), nNeurons);

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)

        fprintf("%.2f, %d \n", cIdx, sIdx)

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

        maxFRs(cIdx, sIdx, :) = max(tuningFns, [], 2);
        
    end
end

save('maxFRs.mat', 'maxFRs')
%% Plot
close all
n = 30;  % pick neuron

[C,S] = meshgrid(contrasts, spreads);
Z = squeeze(maxFRs(:,:,n))';   % match meshgrid dims
figure
surf(C, S, Z)
xlabel('Contrast'); ylabel('Spread'); zlabel('Firing rate')
