clear all
close all
clc

rng('shuffle');

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/ModelV1Uncertainty/DecisionTask/Scripts/')

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
nNeurons     = 200;   % Count of neurons
stimDuration = 2;     % Stimulus duration set to 1 seconds
timeStep     = 0.001; % 0.001s (1ms) - Step size of time bins used for binning stimulus duration 

% Stimulus parameters (angles in radians)
contrasts               = [0.01 0.05]; %[1e-4 0.001 0.01 0.05 0.1 0.2 0.5]; % 0.01
spreads                 = [3 30]; %[3];
uniqStimOris            = 88:0.5:92;
uniqStimOris            = deg2rad(uniqStimOris');
stimParam.numStim       = numel(uniqStimOris);                                % Number of unique stimuli
stimParam.countPerStim  = 60; % 100 
ntrials                 = stimParam.numStim * stimParam.countPerStim * numel(contrasts) * numel(spreads);  % Total number of trials

[c, s, oris] = ndgrid(contrasts, spreads, uniqStimOris);
combinations = [c(:), s(:) oris(:)];
varGain      = getVarGain(combinations(:, 2), combinations(:, 1)); %0.001 * combinations(:, 2) ./ combinations(:, 1);
combinations = [combinations varGain(:)];
trlMatrix    = repmat(combinations, [stimParam.countPerStim 1]);

%% TODO: do this in a loop to make sure its right

% Add random noise to the stimulus orientation
% Stim noise can itself cause gain variability
% Don't forget to shuffle (might not be necessary here though)
trlContrastVector = trlMatrix(:, 1);
trlSpreadVector   = trlMatrix(:, 2);
trlStimVector     = trlMatrix(:, 3);  % Vector of stimuli

gainVector = zeros(size(trlMatrix, 1), nNeurons);
for t = 1:size(trlMatrix, 1)
    gainVector(t,:) = gamrnd(1/trlMatrix(t, 4), trlMatrix(t, 4), [1, nNeurons]);
end

noisyStimVector   = trlStimVector; % Noisy stimulus vector
noisyStimVector   = deg2rad( mod(rad2deg(noisyStimVector), 180) ); % Wrap between 0 and 180

%%
% Shuffle stim vector - this is important for studying effect of top down
% effect
% Shuffle not necessary
% shuffleIdx        = randperm(ntrials);
% trlContrastVector = trlContrastVector(shuffleIdx);
% trlSpreadVector   = trlSpreadVector(shuffleIdx);
% trlStimVector     = trlStimVector(shuffleIdx);
% noisyStimVector   = noisyStimVector(shuffleIdx); 
% gainVector        = gainVector(shuffleIdx, :); % shuffle rows

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
tuningFnSpreads   = linspace(1, 90, 50);

% % TODO: Look at the distribution and verify if the chosen grid spans the
% % 3SD limit (this ensures that the chosen grid is right)
% 
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
%             stimParams.contrastLevel = tuningFnContrasts(cIdx);
%             stimParams.spreadLevel   = deg2rad(tuningFnSpreads(sIdx));
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
%         varGain = getVarGain(tuningFnSpreads(sIdx), tuningFnContrasts(cIdx)); 
%         tuningFnData.data.nrnVarGains.(key1).(key2) = varGain + zeros(1, nNeurons);
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

save('tuningFnData_v2.mat', 'data');

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

for trialIDx = 1:ntrials
    if mod(trialIDx, stimParam.countPerStim/5) == 0
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
    params.contrasts    = tuningFnContrasts; %contrasts;
    params.spreads      = tuningFnSpreads;   %spreads;
    params.uniqStimOris = tuningFnOriSpace;
    
    % Poisson decoder
    [contrastMLE, spreadMLE, thetaMLE, pdfData, MLEs, metrics] = decodePoissonSpikes( ...
        spikes, tuningFns, params);
    decodingError = thetaMLE - trlStimVector(trialIDx); % Error wrt to actual stim ori
    confVar = abs( rad2deg(thetaMLE) - 90) / rad2deg(metrics.sigma);

    trialData.Poisson.pdfData{trialIDx} = pdfData;
    trialData.Poisson.MLEs{trialIDx} = MLEs;
    trialData.Poisson.metrics{trialIDx} = metrics;
    
    % Decoded quantititles
    decodedThetasPossDec(trialIDx)    = thetaMLE;
    decodedContrastsPossDec(trialIDx) = contrastMLE;
    decodedSpreadsPossDec(trialIDx)   = spreadMLE;
    decisionPoissDec(trialIDx)        = rad2deg(thetaMLE) > 90; % CCW if greater than 90
    confVarPoissDec(trialIDx)         = confVar;
    
    % Modulated poisson decoder
    [contrastMLE, spreadMLE, thetaMLE, pdfData, MLEs, metrics] = decodeModulatedPoissonSpikes( ...
        spikes, tuningFns, params, nrnVarGains);
    confVar = abs( rad2deg(thetaMLE) - 90) / rad2deg(metrics.sigma);
    
    trialData.MPoisson.pdfData{trialIDx} = pdfData;
    trialData.MPoisson.MLEs{trialIDx} = MLEs;
    trialData.MPoisson.metrics{trialIDx} = metrics;
    
    % Decoded quantititles
    decodedThetasMPossDec(trialIDx)    = thetaMLE;
    decodedContrastsMPossDec(trialIDx) = contrastMLE;
    decodedSpreadsMPossDec(trialIDx)   = spreadMLE;
    decisionMPoissDec(trialIDx)        = rad2deg(thetaMLE) > 90; % CCW if greater than 90
    confVarMPoissDec(trialIDx)         = confVar;
    
    % Store spike trains and decision results
    %neuronSpikeResponses(trialIDx, :, :) = logical(spikes);  % Store spikes
end

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

%%
data.trialMatrix = trlMatrix;

%%
save('tuningFnData_v2.mat', 'data')

% %%
% % trialData                = tuningFnData.data.trialData;
% % decodedThetasPossDec     = tuningFnData.data.decodedThetasPossDec;
% % decodedContrastsPossDec  = tuningFnData.data.decodedContrastsPossDec;
% % decodedSpreadsPossDec    = tuningFnData.data.decodedSpreadsPossDec;
% % decodedThetasMPossDec    = tuningFnData.data.decodedThetasMPossDec;
% % decodedContrastsMPossDec = tuningFnData.data.decodedContrastsMPossDec;
% % decodedSpreadsMPossDec   = tuningFnData.data.decodedSpreadsMPossDec;
% % neuronSpikeResponses     = tuningFnData.data.neuronSpikeResponses;
% 
% %% Spike - sanity check
% numIntervals = 1;
% intervalSize = floor(length(timeBins) / numIntervals);  % 20 columns per interval
% meanSpkCnt = zeros(nNeurons, numel(contrasts), numel(spreads), numel(uniqStimOris), numIntervals);
% varSpkCnt  = zeros(nNeurons, numel(contrasts), numel(spreads), numel(uniqStimOris), numIntervals);
% 
% for cIdx = 1:numel(contrasts)
%     for sIdx = 1:numel(spreads)
%         for oIdx = 1:numel(uniqStimOris)
%             stimVal  = uniqStimOris(oIdx);
%             trlIdxes = ( trlStimVector ==  stimVal ) & (trlContrastVector == contrasts(cIdx)) & (trlSpreadVector == spreads(sIdx));
%             
%             for nIdx = 1:nNeurons
%                 spkForThisNrn = squeeze( neuronSpikeResponses(trlIdxes, nIdx, :) );
%                 
%                 for i = 1:numIntervals
%                     startCol = (i-1) * intervalSize + 1;
%                     endCol = min(i * intervalSize, length(timeBins));  % Handle the last interval
%                     
%                     % mean spk rate
%                     spkCntInThisInterval = sum( spkForThisNrn(:, startCol:endCol), 2 );
%                     intervalData = mean( spkCntInThisInterval ); 
%                     
%                     meanSpkCnt(nIdx, cIdx, sIdx, oIdx, i) = mean( spkCntInThisInterval );
%                     varSpkCnt(nIdx, cIdx, sIdx, oIdx, i)  = var( spkCntInThisInterval );
%                 end
%             end
%         end
%     end
% end
% 
% % calculate gain values
% nrnVarGainsFit = {};
% 
% for cIdx = 1:numel(contrasts)
%     for sIdx = 1:numel(spreads)
% 
%         gainVals = [];
% 
%         for n=1:nNeurons
%             muSpkCnt     = meanSpkCnt(n, cIdx, sIdx, :);
%             sigma2SpkCnt = varSpkCnt(n, cIdx, sIdx, :);
% 
%             x = muSpkCnt(:); 
%             y = sigma2SpkCnt(:);
% 
%             % Define custom model (edit as needed)
%             ft = fittype('x + sigma_g^2*(x)^2', ...
%                          'independent','x','coefficients',{'sigma_g'});
%             
%             opts = fitoptions(ft);
%             opts.StartPoint = [rand];     % initial guess
%             opts.Lower = [0];          % e.g., sigma_g >= 0
%             % opts.Upper = [10];       % optional upper bound
%             
%             [curve, ~] = fit(x, y, ft, opts);
%             coeffs = coeffvalues(curve);
%             
%             gainVals = [gainVals coeffs];
% 
%         end
% 
%         key1 = sprintf('c_%g', contrasts(cIdx));
%         key2 = sprintf('s_%g', spreads(sIdx));
%         key1 = matlab.lang.makeValidName(key1);
%         key2 = matlab.lang.makeValidName(key2);
%         nrnVarGainsFit.(key1).(key2) = gainVals.^2;
%     end
% end
% 
% figure
% for cIdx = 1:numel(contrasts)
%     for sIdx = 1:numel(spreads)
%         subplot(numel(contrasts), numel(spreads), numel(contrasts)*(cIdx-1) + sIdx )
% 
%         key1 = sprintf('c_%g', contrasts(cIdx));
%         key2 = sprintf('s_%g', spreads(sIdx));
%         key1 = matlab.lang.makeValidName(key1);
%         key2 = matlab.lang.makeValidName(key2);
%         histogram(nrnVarGainsFit.(key1).(key2), DisplayName="Fit")
%         hold on
%         varGain = getVarGain(spreads(sIdx), contrasts(cIdx)); 
%         xline(varGain, DisplayName="GT")
%         hold off
%         title(sprintf("C: %.2f, S: %d", contrasts(cIdx), spreads(sIdx)))
%         xlabel("\sigma_g^2")
%     end
% end
% 
% %% Gain vector for all neuron for specific trials condition
% gainVals = {};
% for trialIDx = 1:ntrials
%     contrastLevel = trlContrastVector(trialIDx); 
%     spreadLevel   = trlSpreadVector(trialIDx);
%     
%     key1 = sprintf('c_%g', contrastLevel);
%     key2 = sprintf('s_%g', spreadLevel);
%     key1 = matlab.lang.makeValidName(key1);
%     key2 = matlab.lang.makeValidName(key2);
% 
%     % Initialize if field doesn't exist
%     if ~isfield(gainVals, key1)
%         gainVals.(key1) = struct();
%     end
%     if ~isfield(gainVals.(key1), key2)
%         gainVals.(key1).(key2) = [];
%     end
% 
%     gainVals.(key1).(key2) = [gainVals.(key1).(key2) gainVector(:,trialIDx)];
% end
% 
% for cIdx = 1:numel(contrasts)
%     for sIdx = 1:numel(spreads)
%         key1 = sprintf('c_%g', contrasts(cIdx));
%         key2 = sprintf('s_%g', spreads(sIdx));
%         key1 = matlab.lang.makeValidName(key1);
%         key2 = matlab.lang.makeValidName(key2);
%         
%         varGain = getVarGain(spreads(sIdx), contrasts(cIdx)); 
% 
%         fprintf("VarGain - Data: %.2f, Actual: %.2f \n", ...
%             var(gainVals.(key1).(key2)(:)), varGain)
%     end
% end
% 
% %%
% figure
% 
% subplot(2, 2, 1)
% scatter(rad2deg(trlStimVector), rad2deg(decodedThetasPossDec), 'filled')
% xlabel("Stim theta (deg)")
% ylabel("Decoded theta (deg)")
% 
% subplot(2, 2, 2)
% scatter(rad2deg(noisyStimVector), rad2deg(decodedThetasPossDec), 'filled')
% xlabel("Stim theta noisy (deg)")
% ylabel("Decoded theta (deg)")
% 
% subplot(2, 2, 3)
% scatter(trlSpreadVector + rand(ntrials, 1), decodedSpreadsPossDec + rand(ntrials, 1), 'filled')
% xlabel("Stim spread (deg)")
% ylabel("Decoded spread (deg)")
% 
% subplot(2, 2, 4)
% scatter(trlContrastVector + 0.001*rand(ntrials, 1), decodedContrastsPossDec + 0.001*rand(ntrials, 1), 'filled')
% xlabel("Stim contrast")
% ylabel("Decoded contrast")
% 
% 
% figure
% 
% subplot(2, 2, 1)
% scatter(rad2deg(trlStimVector), rad2deg(decodedThetasMPossDec), 'filled')
% xlabel("Stim theta (deg)")
% ylabel("Decoded theta (deg)")
% 
% subplot(2, 2, 2)
% scatter(rad2deg(noisyStimVector), rad2deg(decodedThetasMPossDec), 'filled')
% xlabel("Stim theta noisy (deg)")
% ylabel("Decoded theta (deg)")
% 
% subplot(2, 2, 3)
% scatter(trlSpreadVector + rand(ntrials, 1), decodedSpreadsMPossDec + rand(ntrials, 1), 'filled')
% xlabel("Stim spread (deg)")
% ylabel("Decoded spread (deg)")
% 
% subplot(2, 2, 4)
% scatter(trlContrastVector + 0.001*rand(ntrials, 1), decodedContrastsMPossDec + 0.001*rand(ntrials, 1), 'filled')
% xlabel("Stim contrast")
% ylabel("Decoded contrast")
% 
% figure
% scatter(rad2deg(decodedThetasPossDec), rad2deg(decodedThetasMPossDec), 'filled')
% hold on
% plot([0 180], [0 180], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded theta (poisson)")
% ylabel("Decoded theta (modulated poisson)")
% 
% %% Sanity checks
% % 1. mu vs decoded orientation (to ensure that circular mean calculation is correct)
% % 2. Joint vs marginal decoded
% 
% sigmaMPoisson2 = [];
% sigmaPoisson2 = [];
% muPoisson2 = [];
% muMPoisson2 = [];
% 
% sigmaMPoisson = [];
% sigmaPoisson = [];
% muPoisson = [];
% muMPoisson = [];
% 
% thetaMLE_JD_Poisson = []; % Jointly decoded
% thetaMLE_MD_Poisson = []; % Marginally decoded
% thetaMLE_JD_MPoisson = []; % Jointly decoded
% thetaMLE_MD_MPoisson = []; % Marginally decoded
% 
% contrastMLE_JD_Poisson = []; % Jointly decoded
% contrastMLE_MD_Poisson = []; % Marginally decoded
% contrastMLE_JD_MPoisson = []; % Jointly decoded
% contrastMLE_MD_MPoisson = []; % Marginally decoded
% 
% spreadMLE_JD_Poisson = []; % Jointly decoded
% spreadMLE_MD_Poisson = []; % Marginally decoded
% spreadMLE_JD_MPoisson = []; % Jointly decoded
% spreadMLE_MD_MPoisson = []; % Marginally decoded
% 
% for i=1:ntrials
%     sigmaPoisson = [sigmaPoisson tuningFnData.data.trialData.Poisson.metrics{i}.sigma];
%     sigmaMPoisson = [sigmaMPoisson tuningFnData.data.trialData.MPoisson.metrics{i}.sigma];
% 
%     muPoisson = [muPoisson tuningFnData.data.trialData.Poisson.metrics{i}.mu];
%     muMPoisson = [muMPoisson tuningFnData.data.trialData.MPoisson.metrics{i}.mu];
% 
%     sigmaPoisson2 = [sigmaPoisson2 tuningFnData.data.trialData.Poisson.metrics{i}.sigma2];
%     sigmaMPoisson2 = [sigmaMPoisson2 tuningFnData.data.trialData.MPoisson.metrics{i}.sigma2];
%     
%     muPoisson2 = [muPoisson2 tuningFnData.data.trialData.Poisson.metrics{i}.mu2];
%     muMPoisson2 = [muMPoisson2 tuningFnData.data.trialData.MPoisson.metrics{i}.mu2];
%     
%     thetaMLE_JD_Poisson  = [thetaMLE_JD_Poisson tuningFnData.data.trialData.Poisson.MLEs{i}.thetaMLE_v1];
%     thetaMLE_MD_Poisson  = [thetaMLE_MD_Poisson tuningFnData.data.trialData.Poisson.MLEs{i}.thetaMLE_v2];
%     thetaMLE_JD_MPoisson = [thetaMLE_JD_MPoisson tuningFnData.data.trialData.MPoisson.MLEs{i}.thetaMLE_v1];
%     thetaMLE_MD_MPoisson = [thetaMLE_MD_MPoisson tuningFnData.data.trialData.MPoisson.MLEs{i}.thetaMLE_v2];
% 
%     contrastMLE_JD_Poisson  = [contrastMLE_JD_Poisson tuningFnData.data.trialData.Poisson.MLEs{i}.contrastMLE_v1];
%     contrastMLE_MD_Poisson  = [contrastMLE_MD_Poisson tuningFnData.data.trialData.Poisson.MLEs{i}.contrastMLE_v2];
%     contrastMLE_JD_MPoisson = [contrastMLE_JD_MPoisson tuningFnData.data.trialData.MPoisson.MLEs{i}.contrastMLE_v1];
%     contrastMLE_MD_MPoisson = [contrastMLE_MD_MPoisson tuningFnData.data.trialData.MPoisson.MLEs{i}.contrastMLE_v2];
% 
%     spreadMLE_JD_Poisson  = [spreadMLE_JD_Poisson tuningFnData.data.trialData.Poisson.MLEs{i}.spreadMLE_v1];
%     spreadMLE_MD_Poisson  = [spreadMLE_MD_Poisson tuningFnData.data.trialData.Poisson.MLEs{i}.spreadMLE_v2];
%     spreadMLE_JD_MPoisson = [spreadMLE_JD_MPoisson tuningFnData.data.trialData.MPoisson.MLEs{i}.spreadMLE_v1];
%     spreadMLE_MD_MPoisson = [spreadMLE_MD_MPoisson tuningFnData.data.trialData.MPoisson.MLEs{i}.spreadMLE_v2];
% end
% 
% figure
% subplot(2, 2, 1)
% scatter(rad2deg(decodedThetasPossDec), rad2deg(muPoisson), 'filled')
% hold on
% plot([0 180], [0 180], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded theta")
% ylabel("Mu (circular)")
% title("Poisson")
% 
% subplot(2, 2, 2)
% scatter(rad2deg(decodedThetasMPossDec), rad2deg(muMPoisson), 'filled')
% hold on
% plot([0 180], [0 180], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded theta")
% ylabel("Mu (circular)")
% title("Modulated Poisson")
% 
% subplot(2, 2, 3)
% scatter(rad2deg(decodedThetasPossDec), rad2deg(muPoisson2), 'filled')
% hold on
% plot([0 180], [0 180], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded theta")
% ylabel("Mu (linear)")
% title("Poisson")
% 
% subplot(2, 2, 4)
% scatter(rad2deg(decodedThetasMPossDec), rad2deg(muMPoisson2), 'filled')
% hold on
% plot([0 180], [0 180], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded theta")
% ylabel("Mu (linear)")
% title("Modulated Poisson")
% 
% figure
% 
% subplot(2, 2, 1)
% scatter(sigmaPoisson, sigmaMPoisson, 'filled')
% hold on
% plot([0 max(sigmaPoisson)], [0 max(sigmaPoisson)], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Sigma (poisson)")
% ylabel("Sigma (modulated poisson)")
% title("Circular std dev")
% 
% subplot(2, 2, 2)
% scatter(sigmaPoisson2, sigmaMPoisson2, 'filled')
% hold on
% plot([0 max(sigmaPoisson2)], [0 max(sigmaPoisson2)], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Sigma (poisson)")
% ylabel("Sigma (modulated poisson)")
% title("Linear std dev")
% 
% % Joint and marginal decoding
% figure
% 
% subplot(2, 3, 1)
% scatter(rad2deg(thetaMLE_JD_Poisson), rad2deg(thetaMLE_MD_Poisson), 'filled')
% hold on
% plot([0 180], [0 180], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded theta (Joint PDF)")
% ylabel("Decoded theta (Marginal PDF)")
% title("Poisson")
% 
% subplot(2, 3, 4)
% scatter(rad2deg(thetaMLE_JD_MPoisson), rad2deg(thetaMLE_MD_MPoisson), 'filled')
% hold on
% plot([0 180], [0 180], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded theta (Joint PDF)")
% ylabel("Decoded theta (Marginal PDF)")
% title("Modulated Poisson")
% 
% subplot(2, 3, 2)
% scatter((contrastMLE_JD_Poisson), (contrastMLE_MD_Poisson), 'filled')
% hold on
% plot([0 max(contrastMLE_JD_Poisson)], [0 max(contrastMLE_JD_Poisson)], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded contrast (Joint PDF)")
% ylabel("Decoded contrast (Marginal PDF)")
% title("Poisson")
% 
% subplot(2, 3, 5)
% scatter((contrastMLE_JD_MPoisson), (contrastMLE_MD_MPoisson), 'filled')
% hold on
% plot([0 max(contrastMLE_JD_MPoisson)], [0 max(contrastMLE_JD_MPoisson)], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded contrast (Joint PDF)")
% ylabel("Decoded contrast (Marginal PDF)")
% title("Modulated Poisson")
% 
% subplot(2, 3, 3)
% scatter((spreadMLE_JD_Poisson), (spreadMLE_MD_Poisson), 'filled')
% hold on
% plot([0 max(spreadMLE_JD_Poisson)], [0 max(spreadMLE_JD_Poisson)], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded spread (Joint PDF)")
% ylabel("Decoded spread (Marginal PDF)")
% title("Poisson")
% 
% subplot(2, 3, 6)
% scatter((spreadMLE_JD_MPoisson), (spreadMLE_MD_MPoisson), 'filled')
% hold on
% plot([0 max(spreadMLE_JD_MPoisson)], [0 max(spreadMLE_JD_MPoisson)], 'k--', HandleVisibility='off') 
% hold off
% xlabel("Decoded spread (Joint PDF)")
% ylabel("Decoded spread (Marginal PDF)")
% title("Modulated Poisson")
% 
% %% PDFs
% % Randomly pick 25 trial - plot PDF ori for poisson and modulated poisson
% [uniqueRows, ~, groupID] = unique(trlMatrix(:, 1:2), 'rows');
% 
% nGroups = size(uniqueRows,1);
% sampledTrlIdx = cell(nGroups,1);
% 
% for g = 1:nGroups
%     idx = find(groupID == g);           % all trials for this combo
%     sampledTrlIdx{g} = idx(randperm(numel(idx), 5));  % pick 5 random
% end
% 
% figure
% for g = 1:nGroups
%     trlIdxes = sampledTrlIdx{g};
%     
%     for tIdx=1:numel(trlIdxes)
%         
%         varGain = getVarGain(uniqueRows(g, 2), uniqueRows(g, 1));
% 
%         subplot(4, 5, 5*(g-1) + tIdx)
%         hold on
%         plot(rad2deg(tuningFnOriSpace), tuningFnData.data.trialData.Poisson.pdfData{trlIdxes(tIdx)}.pdfOri, ...
%             DisplayName="Poisson", LineWidth=1.5)
%         plot(rad2deg(tuningFnOriSpace), tuningFnData.data.trialData.MPoisson.pdfData{trlIdxes(tIdx)}.pdfOri, ...
%             DisplayName="Modulated Poisson", LineWidth=1.5)
%         %xline(tuningFnData.data.trialData.Poisson.metrics{trlIdxes(tIdx)}.mu, LineWidth=1.5, LineStyle="--")
%         hold off
%         xlabel("orientation (deg)")
%         ylabel("P(ori)")
%         legend
%         xStart = rad2deg(tuningFnData.data.trialData.Poisson.metrics{trlIdxes(tIdx)}.mu) - 20;
%         xEnd   = rad2deg(tuningFnData.data.trialData.Poisson.metrics{trlIdxes(tIdx)}.mu) + 20;
%         xlim([xStart, xEnd])
%         title(sprintf("C: %.2f, S: %d, varGain: %.2f", ...
%             uniqueRows(g, 1), uniqueRows(g, 2), varGain))
% 
%     end
% end
% 
% %%
% % function varGain = getVarGain(spread, contrast)
% %     varGain = 0.0001 * spread ./ contrast;
% % end
% 
% % tuningFnContrasts = linspace(1e-4, 0.15, 49);
% % tuningFnSpreads   = linspace(1, 90, 50);
% % % tuningFnContrasts = linspace(0.01, 0.05, 10);
% % % tuningFnSpreads   = linspace(3, 30, 10);
% % 
% % gains = zeros(numel(tuningFnContrasts), numel(tuningFnSpreads));
% % 
% % for i=1:numel(tuningFnContrasts)
% %     for j=1:numel(tuningFnSpreads)
% %         gains(i, j) = getVarGain(tuningFnSpreads(j), tuningFnContrasts(i));
% %     end
% % end
% % 
% % figure
% % plot(tuningFnSpreads, gains(1, :))
% % hold on
% % plot(tuningFnSpreads, gains(end, :))
% % hold off
% % % ylim([0 1])

function varGain = getVarGain(spread, contrast)
%     sigmaG = (0.15/27)*spread + (0.25 - (0.15/27)*30 ) + ...
%         ( - (0.2/0.04)*contrast + (0.3 + (0.2/0.04)*0.01) + 0.4);
    sigmaG = (0.15/27)*spread + (0.25 - (0.15/27)*30 ) + ...
        ( - (0.2/0.04)*contrast + (0.3 + (0.2/0.04)*0.01) + 0.4);
    varGain = sigmaG.^2;
end

% % baselines 0.0225, 0.0625
% % peaks: 0.2025, 0.0784