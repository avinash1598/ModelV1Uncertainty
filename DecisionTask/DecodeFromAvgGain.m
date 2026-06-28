clear all
close all
clc

rng('shuffle');
warning('Make sure spread used is in rad')
addpath('C:\Users\avinash1598\Desktop\ModelV1Uncertainty\DecisionTask\Scripts\')

% tuningFnData = load('tuningFnDataGainMatrix.mat');
tuningFnData = load('tuningFnDataDecodeAvgGain.mat');

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
contrasts               = [0.01 0.05]; %[1e-4 0.001 0.01 0.05 0.1 0.2 0.5]; % 0.01
spreads                 = deg2rad([3 30]); %[3];
uniqStimOris            = linspace(90-10, 90+10, 21); %80:0.5:100;
uniqStimOris            = deg2rad(uniqStimOris');
stimParam.numStim       = numel(uniqStimOris);                                % Number of unique stimuli
stimParam.countPerStim  = 60; % 100 
ntrials                 = stimParam.numStim * stimParam.countPerStim * numel(contrasts) * numel(spreads);  % Total number of trials

[c, s, oris] = ndgrid(contrasts, spreads, uniqStimOris);
combinations = [c(:), s(:) oris(:)];
trlMatrix    = repmat(combinations, [stimParam.countPerStim 1]);

%% TODO: do this in a loop to make sure its right

% Add random noise to the stimulus orientation
% Stim noise can itself cause gain variability
% Don't forget to shuffle (might not be necessary here though)
trlContrastVector = trlMatrix(:, 1);
trlSpreadVector   = trlMatrix(:, 2);
trlStimVector     = trlMatrix(:, 3);  % Vector of stimuli

gainVector = zeros(size(trlMatrix, 1), nNeurons);
varGains   = zeros(size(trlMatrix, 1), nNeurons);
for t = 1:size(trlMatrix, 1)
    varGain = getVarGainForAllNeuron(nNeurons, ...
        trlSpreadVector(t), trlContrastVector(t));
    gainVector(t,:) = gamrnd(1./varGain, varGain);
    varGains(t, :)  = varGain;
end

avgGainVector = zeros(size(trlMatrix, 1), nNeurons);
meanVarGains   = zeros(size(trlMatrix, 1), nNeurons);
for t = 1:size(trlMatrix, 1)
    varGainMean = getVarGain(trlSpreadVector(t), trlContrastVector(t)); 
    avgGainVector(t,:) = gamrnd(1/varGainMean, varGainMean, [1, nNeurons]);
    meanVarGains(t, :) = varGainMean;
end

stimNoise         = 0 + 0.1 * randn(1, ntrials); % What is std dev here? 5.73 degrees
noisyStimVector   = trlStimVector + stimNoise; % Noisy stimulus vector
noisyStimVector   = deg2rad( mod(rad2deg(noisyStimVector), 180) ); % Wrap between 0 and 180

figure
idx = 1;
for i = [10, 237, 1073, 2000]
    subplot(3, 4, idx)
    histogram(sqrt( varGains(i,:) ), BinWidth=0.05, DisplayName="Variable \sigma_g")
    hold on
    xline(mean( sqrt( meanVarGains(i,:) ) ), lineStyle="--", lineWidth=2, DisplayName="Avg \sigma_g")
    hold off
    xlabel('\sigma_g')
    ylabel('Count')
    title('\sigma_g' + sprintf(" distribution (across neurons) \nStimulus Difficulty level %d", idx))
    xlim([0 2])
    legend
    
    subplot(3, 4, 4 + idx)
    histogram(gainVector(i,:), BinWidth=0.5, DisplayName="Variable \sigma_g")
    hold on
    histogram(avgGainVector(i,:), BinWidth=0.5, DisplayName="Avg \sigma_g")
    hold off
    xlabel('Gains')
    ylabel('Count')
    title('Gain samples')
    xlim([0 6])
    legend

    subplot(3, 4, 8 + idx)
    scatter(sqrt(varGains(i, :)), sqrt(varGains(i+1, :)))
    xlabel('\sigma_g (Difficulty level 1)')
    ylabel('\sigma_g (Difficulty level 2)')
    title('Covarying gains for neurons')
    legend

    idx = idx + 1;
end

%%
% % Neurons tuning parameters
% tuningParams.d     = zeros(nNeurons, 1) + 0;   % (fixed) Direction selectivity - set it to zero (no need for neuron to be directional selective).
% tuningParams.alpha = zeros(nNeurons, 1) + 2;   % (fixed) Aspect ratio - controls sharpness. Keep this fixed. Reducing the value makes the changes very rapid towards the end which we probably don't want.
% tuningParams.b     = zeros(nNeurons, 1) + 2;   % (fixed maybe/variable - (0.5, some max - 3, 4 ...)) Control this - Control sharpness + range of the neuron. Set it to 2 for these simulations
% tuningParams.q     = zeros(nNeurons, 1) + 2;   % 1 or 2?? (variable) Set it to some constant. Controls the sharpness and amplitude of peak FR.
% tuningParams.w     = zeros(nNeurons, 1) + 0;   % 0 or 1?? Doesn't matter since untuned component is zero. Weight of untuned filter. Set it to 0. (fixed) Doesn't matter what is val is becz untuned filter amp is zero.
% tuningParams.e1    = zeros(nNeurons, 1);       % Stimulus independent spontaneous discharge (variable) Controls dynamic range.
% tuningParams.e2    = zeros(nNeurons, 1);       % Stimulus dependent spontaneous discharge (variable) Controls dynamic range.
% tuningParams.gam   = zeros(nNeurons, 1);       % Controls response amplitude (variable) Controls dynamic range.
% tuningParams.beta  = zeros(nNeurons, 1) + 0;   % Stimulus independent constant (variable) Controls dynamic range.
% tuningParams.UNTUNED_FILTER_AMPL = 0;          % (fixed) Untuned filter not needed.
% tuningParams.eps   = 0;                        % (this probably needs to be sampled every trial??) Not right place to update here Normalization noise sampled from some distribution with sigma_g standard deviation
% 
% tuningParams.alpha(:)   = 2; % 1 + 5*rand(1, nNeurons);  Aspect ratio uniformly sampled from 0 - 5
% tuningParams.b(:)       = 2; % Derivative order - Just set it 2 (Zoey's paper). Non integer vaalues can give imaginary values. More like log uniform 0.0125 + 8*rand(1, nNeurons); 
% tuningParams.q(:)       = exp( log(1.8)*rand(nNeurons, 1) );  % Transduction (this is not uniformly distributed) exp( log(1.8)*rand(1, nNeurons) )
% tuningParams.beta(:)    = lognrnd(2.5, 0.5, nNeurons, 1);     % lognrnd(2.5, 0.5, 1, nNeurons);     % Normalization constant lognrnd(mu, sigma, 1, nneurons)
% 
% tuningParams.e1(:)        = lognrnd(0.8, 0.6, nNeurons, 1);       % Range: 0 - 10 ips, Might not be correct initilization
% tuningParams.gam(:)       = 4000;
% 
% neuronsPrefOrientation    = zeros(nNeurons, 1);
% neuronsPrefOrientation(:) = pi * rand(nNeurons, 1);        % Randomly choose neurons preferred orientation from 0 to pi (not directional selective)

timeBins                  = 0:timeStep:stimDuration; 
stimRespProfile           = 1 + zeros(1, numel(timeBins));
gainVector                = gainVector'; %nNeurons x nTrials 1 + zeros(nNeurons, ntrials); % constant gain - NO gain modulation

tuningParams              = tuningFnData.data.tuningParams;
neuronsPrefOrientation    = tuningFnData.data.neuronsPrefOrientation;

% tuningFns                 = tuningFnData.data.tuningFns;
% nrnVarGains               = tuningFnData.data.nrnVarGains;

%% Update tuning function data
tuningFnOriSpace  = linspace(0, pi, 361);
tuningFnContrasts = linspace(1e-4, 0.15, 49);
tuningFnSpreads   = deg2rad(linspace(1, 90, 50));

% % % TODO: Look at the distribution and verify if the chosen grid spans the
% % % 3SD limit (this ensures that the chosen grid is right)
% for cIdx = 1:numel(tuningFnContrasts)
%     for sIdx = 1:numel(tuningFnSpreads)
% 
%         fprintf("%.2f, %d \n", cIdx, sIdx)
% 
%         tuningFns = [];
% 
%         for oIdx = 1:numel(tuningFnOriSpace)
% 
%             % Get tuning params every trial
%             stimParams = {};
%             stimParams.contrastLevel = tuningFnContrasts(cIdx);
%             stimParams.spreadLevel   = tuningFnSpreads(sIdx);
%             stimParams.stimOri       = tuningFnOriSpace(oIdx);
% 
%             tFn = getOriTunedStimRespFunction( ...
%                     neuronsPrefOrientation, tuningParams, stimParams);
%             firingRates = tFn.FR; % Firing rate for this trial
% 
%             tuningFns = [tuningFns firingRates];
%         end
% 
%         key1 = sprintf('c_%g', tuningFnContrasts(cIdx));
%         key2 = sprintf('s_%g', tuningFnSpreads(sIdx));
%         key1 = matlab.lang.makeValidName(key1);
%         key2 = matlab.lang.makeValidName(key2);
%         tuningFnData.data.tuningFns.(key1).(key2) = tuningFns;
% 
%         % gains
%         %varGain = getVarGainForAllNeuron(nNeurons, ...
%         %    tuningFnSpreads(sIdx), tuningFnContrasts(cIdx));
%         % This will be used for decoding. Use average sigma_g 
%         % tuningFnData.data.nrnVarGains.(key1).(key2) = varGain'; % make sure the size is 1xnNeurons
%         varGain = getVarGain(tuningFnSpreads(sIdx), tuningFnContrasts(cIdx)); 
%         tuningFnData.data.nrnVarGains.(key1).(key2) = varGain + zeros(1, nNeurons);
% 
%     end
% end

data.stimDuration = stimDuration;    
data.timeStep     = timeStep; 
data.nNeurons     = nNeurons;
data.contrasts    = contrasts;
data.spreads      = spreads;
data.uniqStimOris = uniqStimOris;
data.tuningParams = tuningParams;
data.neuronsPrefOrientation = neuronsPrefOrientation;
data.tuningFns    = tuningFnData.data.tuningFns;
data.nrnVarGains  = tuningFnData.data.nrnVarGains;

% Update tuning Fns
tuningFns    = tuningFnData.data.tuningFns;
nrnVarGains  = tuningFnData.data.nrnVarGains;

% save('tuningFnDataDecodeAvgGain.mat', 'data');

%% Generate spike data
% Structures to store final neuron spikes
% Preallocate result and spike response matrices
% firingRates          = squeeze( tuningFn(:, 1, 1, :) );
decodedThetasPossDec        = zeros(ntrials, 1);
decodedContrastsPossDec     = zeros(ntrials, 1);
decodedSpreadsPossDec       = zeros(ntrials, 1);
decisionPoissDec            = zeros(ntrials, 1);
confVarPoissDec             = zeros(ntrials, 1);

decodedThetasMPossDec        = zeros(ntrials, 1);
decodedContrastsMPossDec     = zeros(ntrials, 1);
decodedSpreadsMPossDec       = zeros(ntrials, 1);
decisionMPoissDec            = zeros(ntrials, 1);
confVarMPoissDec             = zeros(ntrials, 1);

%neuronSpikeResponses = false(ntrials, nNeurons, length(timeBins)); 

trialData = {};

% Temporary data structures
poissPdf     = cell(ntrials,1);
poissMLEs    = cell(ntrials,1);
poissMetrics = cell(ntrials,1);

mPoissPdf     = cell(ntrials,1);
mPoissMLEs    = cell(ntrials,1);
mPoissMetrics = cell(ntrials,1);

p = gcp('nocreate');   
if isempty(p)
    parpool;           
end

parfor trialIDx = 1:ntrials
    if mod(trialIDx, 60) == 0
        disp(trialIDx)
    end
    
    % Get tuning params every trial
    stimParams = {};
    stimParams.contrastLevel = trlContrastVector(trialIDx); %0.01 0.04 0.1 Use these two values. Saturation happens pretty quickly
    %stimParams.spreadLevel   = deg2rad(trlSpreadVector(trialIDx));
    stimParams.spreadLevel   = trlSpreadVector(trialIDx);
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
    params.contrasts    = tuningFnContrasts; %contrasts;
    params.spreads      = tuningFnSpreads;   %spreads;
    params.uniqStimOris = tuningFnOriSpace;

    % Poisson decoder
    [contrastMLE, spreadMLE, thetaMLE, pdfData, MLEs, metrics] = decodePoissonSpikes( ...
        spikes, tuningFns, params);
    decodingError = thetaMLE - trlStimVector(trialIDx); % Error wrt to actual stim ori
    confVar = abs( rad2deg(thetaMLE) - 90) / rad2deg(metrics.sigma);

    poissPdf{trialIDx} = pdfData;
    poissMLEs{trialIDx} = MLEs;
    poissMetrics{trialIDx} = metrics;
    
    % Decoded quantititles
    decodedThetasPossDec(trialIDx)    = thetaMLE;
    decodedContrastsPossDec(trialIDx) = contrastMLE;
    decodedSpreadsPossDec(trialIDx)   = spreadMLE;
    decisionPoissDec(trialIDx)        = rad2deg(thetaMLE) > 90; % CCW if greater than 90
    confVarPoissDec(trialIDx)         = confVar;
    
    %keyboard

    % Modulated poisson decoder
    [contrastMLE, spreadMLE, thetaMLE, pdfData, MLEs, metrics] = decodeModulatedPoissonSpikes( ...
        spikes, tuningFns, params, nrnVarGains);
    confVar = abs( rad2deg(thetaMLE) - 90) / rad2deg(metrics.sigma);
    
    mPoissPdf{trialIDx} = pdfData;
    mPoissMLEs{trialIDx} = MLEs;
    mPoissMetrics{trialIDx} = metrics;
    
    % Decoded quantititles
    decodedThetasMPossDec(trialIDx)    = thetaMLE;
    decodedContrastsMPossDec(trialIDx) = contrastMLE;
    decodedSpreadsMPossDec(trialIDx)   = spreadMLE;
    decisionMPoissDec(trialIDx)        = rad2deg(thetaMLE) > 90; % CCW if greater than 90
    confVarMPoissDec(trialIDx)         = confVar;
    
    % Store spike trains and decision results
    %neuronSpikeResponses(trialIDx, :, :) = logical(spikes);  % Store spikes
end

trialData.Poisson.pdfData = poissPdf;
trialData.Poisson.MLEs = poissMLEs;
trialData.Poisson.metrics = poissMetrics;

trialData.MPoisson.pdfData = mPoissPdf;
trialData.MPoisson.MLEs = mPoissMLEs;
trialData.MPoisson.metrics = mPoissMetrics;

% tuningFnData.data.trialData = trialData;
% tuningFnData.data.decodedThetasPossDec = decodedThetasPossDec;
% tuningFnData.data.decodedContrastsPossDec = decodedContrastsPossDec;
% tuningFnData.data.decodedSpreadsPossDec = decodedSpreadsPossDec;
% tuningFnData.data.decodedThetasMPossDec = decodedThetasMPossDec;
% tuningFnData.data.decodedContrastsMPossDec = decodedContrastsMPossDec;
% tuningFnData.data.decodedSpreadsMPossDec = decodedSpreadsMPossDec;
% tuningFnData.data.neuronSpikeResponses = neuronSpikeResponses;

%%
data.trialData                = trialData;
data.decodedThetasPossDec     = decodedThetasPossDec;
data.decodedContrastsPossDec  = decodedContrastsPossDec;
data.decodedSpreadsPossDec    = decodedSpreadsPossDec;
data.decodedThetasMPossDec    = decodedThetasMPossDec;
data.decodedContrastsMPossDec = decodedContrastsMPossDec;
data.decodedSpreadsMPossDec   = decodedSpreadsMPossDec;
%data.neuronSpikeResponses    = neuronSpikeResponses;
data.decisionPoissDec         = decisionPoissDec;
data.confVarPoissDec          = confVarPoissDec;
data.decisionMPoissDec        = decisionMPoissDec;
data.confVarMPoissDec         = confVarMPoissDec;

data.trialMatrix = trlMatrix;
data.gainVector  = gainVector;
data.noisyStimVector = noisyStimVector;

%%
save('tuningFnDataDecodeAvgGain.mat', 'data')

disp("Done")
%%
function nrnsVarGain = getVarGainForAllNeuron(nNeurons, spread, contrast)
    sigmaG = sqrt(getVarGain(spread, contrast));

    m = sigmaG; % Note relationship is wrt sigma_g (not sigma_g2)
    v = ( 0.4*m ).^2;
    sigma2 = log(1 + v./m.^2);
    sigma = sqrt(sigma2);
    mu = log(m) - 0.5*sigma2;
    
    X = lognrnd(mu, sigma, nNeurons, 1);
    nrnsVarGain = X.^2;
    
    % Should be pretty close
    assert(abs(sigmaG - mean(sqrt(nrnsVarGain))) < 0.1)
end

function varGain = getVarGain(spread, contrast)
    spread = rad2deg(spread);
    sigmaG = (0.15/27)*spread + (0.25 - (0.15/27)*30 ) + ...
        ( - (0.2/0.04)*contrast + (0.3 + (0.2/0.04)*0.01) + 0.4);
    varGain = sigmaG.^2;
end