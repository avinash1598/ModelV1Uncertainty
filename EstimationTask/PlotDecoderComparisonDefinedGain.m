clear all
close all
clc

rng('shuffle');

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/ModelV1Uncertainty/EstimationTask/Scripts/')

set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesLineWidth',1, ...
    'defaultAxesBox','off', ...
    'defaultAxesTickDir','out', ...
    'defaultLineLineWidth',1.5, ...
    'defaultAxesLabelFontSizeMultiplier',1.1, ...
    'defaultAxesTitleFontWeight','normal', ...
    'defaultLegendBox','off');

% tuningFnData = load('tuningFnData.mat');
% tuningFnData = load('tuningFnData_v2.mat');
tuningFnData = load('tuningFnData_v3.mat');

% ----------------------------------
% Params
% ----------------------------------
nNeurons     = tuningFnData.data.nNeurons;      % Count of neurons
stimDuration = tuningFnData.data.stimDuration;  % Stimulus duration set to 1 seconds
timeStep     = tuningFnData.data.timeStep;      % 0.001s (1ms) - Step size of time bins used for binning stimulus duration 

% Stimulus parameters (angles in radians)
contrasts               = tuningFnData.data.contrasts; 
spreads                 = tuningFnData.data.spreads;
uniqStimOris            = tuningFnData.data.uniqStimOris;
% uniqStimOris            = deg2rad(uniqStimOris');
stimParam.numStim       = numel(uniqStimOris);                                
stimParam.countPerStim  = 20; % Don't change this 
ntrials                 = stimParam.numStim * stimParam.countPerStim * numel(contrasts) * numel(spreads);  

[c, s, oris] = ndgrid(contrasts, spreads, uniqStimOris);
combinations = [c(:), s(:) oris(:)];
varGain      = getVarGain(combinations(:, 2), combinations(:, 1)); %0.001 * combinations(:, 2) ./ combinations(:, 1);
combinations = [combinations varGain(:)];
trlMatrix    = repmat(combinations, [stimParam.countPerStim 1]);

%%
trlContrastVector = trlMatrix(:, 1);
trlSpreadVector   = trlMatrix(:, 2);
trlStimVector     = trlMatrix(:, 3);  % Vector of stimuli

gainVector = zeros(size(trlMatrix, 1), nNeurons);
for t = 1:size(trlMatrix, 1)
    gainVector(t,:) = gamrnd(1/trlMatrix(t, 4), trlMatrix(t, 4), [1, nNeurons]);
end

noisyStimVector   = trlStimVector; % Noisy stimulus vector
noisyStimVector   = deg2rad( mod(rad2deg(noisyStimVector), 180) ); % Wrap between 0 and 180

timeBins                  = 0:timeStep:stimDuration; 
stimRespProfile           = 1 + zeros(1, numel(timeBins));
gainVector                = gainVector'; %nNeurons x nTrials 1 + zeros(nNeurons, ntrials); % constant gain - NO gain modulation
tuningParams              = tuningFnData.data.tuningParams;
neuronsPrefOrientation    = tuningFnData.data.neuronsPrefOrientation;
tuningFns                 = tuningFnData.data.tuningFns;
nrnVarGains               = tuningFnData.data.nrnVarGains;

% Don't change it
tuningFnOriSpace          = linspace(0, pi, 361);
% tuningFnContrasts         = linspace(1e-4, 0.2, 49); % v1
tuningFnContrasts         = linspace(1e-4, 0.15, 49); % v2
tuningFnSpreads           = linspace(1, 90, 50); 

%% Tuning fns for given stimulus condition
tuningFnsByExpConditions = {};

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)

        tuningFns_ = [];
        
        for oIdx = 1:numel(tuningFnOriSpace)
            
            % Get tuning params every trial
            stimParams.contrastLevel = contrasts(cIdx);
            stimParams.spreadLevel   = deg2rad(spreads(sIdx));
            stimParams.stimOri       = tuningFnOriSpace(oIdx);
            
            tFn = getOriTunedStimRespFunction( ...
                    neuronsPrefOrientation, tuningParams, stimParams);
            firingRates = tFn.FR; % Firing rate for this trial
            
            tuningFns_ = [tuningFns_ firingRates];
        end

        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);
        tuningFnsByExpConditions.(key1).(key2) = tuningFns_;
    end
end

%% Spike related data
trialData                = tuningFnData.data.trialData;
decodedThetasPossDec     = tuningFnData.data.decodedThetasPossDec;
decodedContrastsPossDec  = tuningFnData.data.decodedContrastsPossDec;
decodedSpreadsPossDec    = tuningFnData.data.decodedSpreadsPossDec;
decodedThetasMPossDec    = tuningFnData.data.decodedThetasMPossDec;
decodedContrastsMPossDec = tuningFnData.data.decodedContrastsMPossDec;
decodedSpreadsMPossDec   = tuningFnData.data.decodedSpreadsMPossDec;
neuronSpikeResponses     = tuningFnData.data.neuronSpikeResponses;

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

%% data by stim condition
pdfOriByCondition = {};
estimationErr = {};
sigma_theta = {};
decodedThetaDiff = {};
decodedSigmaDiff = {};

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        estimationErr.MPoisson.(key1).(key2) = [];
        estimationErr.Poisson.(key1).(key2)  = [];
        sigma_theta.MPoisson.(key1).(key2)   = [];
        sigma_theta.Poisson.(key1).(key2)    = [];
        decodedThetaDiff.(key1).(key2)       = [];
        decodedSigmaDiff.(key1).(key2)       = [];
    end
end

for i=1:ntrials
    key1 = sprintf('c_%g', trlMatrix(i, 1));
    key2 = sprintf('s_%g', trlMatrix(i, 2));
    key3 = sprintf('o_%g', trlMatrix(i, 3));
    key1 = matlab.lang.makeValidName(key1);
    key2 = matlab.lang.makeValidName(key2);
    key3 = matlab.lang.makeValidName(key3);
    
    pdfOri = tuningFnData.data.trialData.MPoisson.pdfData{i}.pdfOri;
    pdfOriByCondition.MPoisson.(key1).(key2).(key3) = pdfOri;

    pdfOri = tuningFnData.data.trialData.Poisson.pdfData{i}.pdfOri;
    pdfOriByCondition.Poisson.(key1).(key2).(key3) = pdfOri;

    % Modulated Poisson
    err = rad2deg(decodedThetasMPossDec(i)) - rad2deg(noisyStimVector(i));
    err = mod(err + 90, 180) - 90; % Error between -90 and 90
    estimationErr.MPoisson.(key1).(key2) = [estimationErr.MPoisson.(key1).(key2) err];
    sigma_theta.MPoisson.(key1).(key2) = [sigma_theta.MPoisson.(key1).(key2) ...
        trialData.MPoisson.metrics{i}.sigma];

    % Poisson
    err = rad2deg(decodedThetasPossDec(i)) - rad2deg(noisyStimVector(i));
    err = mod(err + 90, 180) - 90; % Error between -90 and 90
    estimationErr.Poisson.(key1).(key2) = [estimationErr.Poisson.(key1).(key2) err];
    sigma_theta.Poisson.(key1).(key2) = [sigma_theta.Poisson.(key1).(key2) ...
        trialData.Poisson.metrics{i}.sigma];

    decodedThetaDiff.(key1).(key2) = [decodedThetaDiff.(key1).(key2) ... 
        rad2deg( (decodedThetasMPossDec(i) - decodedThetasPossDec(i)) )];

    decodedSigmaDiff.(key1).(key2) = [decodedSigmaDiff.(key1).(key2) ...
        rad2deg( ( tuningFnData.data.trialData.MPoisson.metrics{i}.sigma - ...
    tuningFnData.data.trialData.Poisson.metrics{i}.sigma))];
end

%% Confidence for each trial
% 1. Assuming modulated poisson as ground truth of behavior
% 2. Assuming poisson as ground truth of behavior
% Assumption: choose confidence criteria which gives 50:50 HC and LC split

sigmaPoisson  = [];
sigmaMPoisson = [];

for i=1:ntrials
    % Note: 17 values are NaN: Not a big problem but keep in mind
    sigmaPoisson  = [sigmaPoisson tuningFnData.data.trialData.Poisson.metrics{i}.sigma];
    sigmaMPoisson = [sigmaMPoisson tuningFnData.data.trialData.MPoisson.metrics{i}.sigma];
end

trlFltIdx = sigmaMPoisson < 1; %0.04;
sigmaMPoisson = sigmaMPoisson(trlFltIdx);
sigmaPoisson  = sigmaPoisson(trlFltIdx);

p = 50;
Cc_MPoiss = prctile(sigmaMPoisson, p);
Cc_Poiss  = prctile(sigmaPoisson,  p);
% Cc_MPoiss = median(sigmaMPoisson); % greater than this low confidence, less than this high confidence
% Cc_Poiss  = median(sigmaPoisson);  % greater than this low confidence, less than this high confidence

confReports_MPoiss = sigmaMPoisson < Cc_MPoiss; 
confReports_Poiss  = sigmaPoisson < Cc_Poiss;

%% confidecence report vs sigma values - logistic regression
% For recovery it's important that there is good split of high confidence
% and low confidence trials. All HC or all LC implies no behavioral
% variability and hence recovery becomes impossible. 

warning('off','all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ground truth: Modulated Poisson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bA = glmfit(sigmaPoisson(:), confReports_MPoiss(:), 'binomial');
pA = glmval(bA, sigmaPoisson(:), 'logit'); % todo maybe provide test data, P(confReports=1∣σ_MPoisson​=x)

bB = glmfit(sigmaMPoisson(:), confReports_MPoiss(:), 'binomial');
pB = glmval(bB, sigmaMPoisson(:), 'logit'); % Here probability is for high confidence given the uncertainty estimate

% Plotting business
figure;

xFit = linspace(min(sigmaPoisson), max(sigmaPoisson), 200);
pFit = glmval(bA, xFit, 'logit');

subplot(2, 2, 1)
scatter(sigmaPoisson, confReports_MPoiss, 20, 'filled', DisplayName='data'); hold on;
plot(xFit, pFit, 'LineWidth', 2, DisplayName='Logistic Fit');
xlabel('\sigma_{Poisson}');
ylabel('P(Conf = 1)');
ylim([-0.05 1.05]);
xlim([0 0.18])
legend
title("GT (conf labels): M Poisson")


xFit = linspace(min(sigmaMPoisson), max(sigmaMPoisson), 200);
pFit = glmval(bB, xFit, 'logit');

subplot(2, 2, 3)
scatter(sigmaMPoisson, confReports_MPoiss, 20, 'filled'); hold on;
plot(xFit, pFit, 'LineWidth', 2);
xlabel('\sigma_{Modulated Poisson}');
ylabel('P(Conf = 1)');
ylim([-0.05 1.05]);
xlim([0 0.18])
title("GT (conf labels): M Poisson")


% Likelihood P(conf report/sigma_theta) = 
% for each trial likelihood = P(confReports=1∣σ_MPoisson​=x) if conf == 1, 
% = 1 - P(confReports=1∣σ_MPoisson​=x) if conf == 0
LL_A = sum(confReports_MPoiss(:).*log(pA(:)) + (1-confReports_MPoiss(:)).*log(1-pA(:)));
LL_B = sum(confReports_MPoiss(:).*log(pB(:)) + (1-confReports_MPoiss(:)).*log(1-pB(:)));

% % Log-likelihood
% pA = glmval(bA, U1, 'logit');
% LL_A = sum(y.*log(pA) + (1-y).*log(1-pA));

recoveryLL_GT_MPoiss = [LL_A LL_B]; % Poiss, MPoiss

fprintf("Log Likelihood comparison\n\n")
fprintf("----------------------------------------\n " );
fprintf("Ground truth: Modulated Poisson decoder\n" );
fprintf("----------------------------------------\n" );
fprintf("Poisson decoder: %.4f \n", LL_A );
fprintf("Modulated Poisson decoder: %.4f \n", LL_B );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ground truth: Poisson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bA = glmfit(sigmaPoisson(:), confReports_Poiss(:), 'binomial');
pA = glmval(bA, sigmaPoisson(:), 'logit'); % todo maybe provide test data, P(confReports=1∣σ_MPoisson​=x)

bB = glmfit(sigmaMPoisson(:), confReports_Poiss(:), 'binomial');
pB = glmval(bB, sigmaMPoisson(:), 'logit'); % Here probability is for high confidence given the uncertainty estimate

% Likelihood P(conf report/sigma_theta) = 
% for each trial likelihood = P(confReports=1∣σ_MPoisson​=x) if conf == 1, 
% = 1 - P(confReports=1∣σ_MPoisson​=x) if conf == 0
LL_A = sum(confReports_Poiss(:).*log(pA(:)) + (1-confReports_Poiss(:)).*log(1-pA(:)));
LL_B = sum(confReports_Poiss(:).*log(pB(:)) + (1-confReports_Poiss(:)).*log(1-pB(:)));

% % Log-likelihood
% pA = glmval(bA, U1, 'logit');
% LL_A = sum(y.*log(pA) + (1-y).*log(1-pA));

fprintf("Log Likelihood comparison\n\n")
fprintf("----------------------------------------\n " );
fprintf("Ground truth: Poisson decoder\n" );
fprintf("----------------------------------------\n" );
fprintf("Poisson decoder: %.4f \n", LL_A );
fprintf("Modulated Poisson decoder: %.4f \n", LL_B );

recoveryLL_GT_Poiss = [LL_A LL_B]; % Poiss, MPoiss

xFit = linspace(min(sigmaPoisson), max(sigmaPoisson), 200);
pFit = glmval(bA, xFit, 'logit');

subplot(2, 2, 2)
scatter(sigmaPoisson, confReports_Poiss, 20, 'filled', DisplayName='data'); hold on;
plot(xFit, pFit, 'LineWidth', 2, DisplayName='Logistic Fit');
xlabel('\sigma_{Poisson}');
ylabel('P(Conf = 1)');
ylim([-0.05 1.05]);
xlim([0 0.18])
legend
title("GT (conf labels): Poisson")

xFit = linspace(min(sigmaMPoisson), max(sigmaMPoisson), 200);
pFit = glmval(bB, xFit, 'logit');

subplot(2, 2, 4)
scatter(sigmaMPoisson, confReports_Poiss, 20, 'filled'); hold on;
plot(xFit, pFit, 'LineWidth', 2);
xlabel('\sigma_{Modulated Poisson}');
ylabel('P(Conf = 1)');
ylim([-0.05 1.05]);
xlim([0 0.18])
title("GT (conf labels): Poisson")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross-validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ground truth: Modulated Poisson
k = 5;
cv = cvpartition(confReports_MPoiss(:), 'KFold', k);

lossA = crossval(@(Xtrain,Ytrain,Xtest,Ytest) ...
    logLoss(glmfit(Xtrain, Ytrain, 'binomial'), Xtest, Ytest), ...
    sigmaPoisson(:), confReports_MPoiss(:), 'partition', cv);

lossB = crossval(@(Xtrain,Ytrain,Xtest,Ytest) ...
    logLoss(glmfit(Xtrain, Ytrain, 'binomial'), Xtest, Ytest), ...
    sigmaMPoisson(:), confReports_MPoiss(:), 'partition', cv);

fprintf("\n\nCross-validation likelihood comparison\n\n")
fprintf("----------------------------------------\n " );
fprintf("Ground truth: Modulated Poisson decoder \n" );
fprintf("----------------------------------------\n" );
fprintf("Poisson decoder: %.4f \n", -mean(lossA) );
fprintf("Modulated Poisson decoder: %.4f \n", -mean(lossB) );

cvLL_GT_MPoiss = [-lossA -lossB]; % Poiss, MPoiss

% Ground truth: Poisson
k = 5;
cv = cvpartition(confReports_Poiss(:), 'KFold', k);

lossA = crossval(@(Xtrain,Ytrain,Xtest,Ytest) ...
    logLoss(glmfit(Xtrain, Ytrain, 'binomial'), Xtest, Ytest), ...
    sigmaPoisson(:), confReports_Poiss(:), 'partition', cv);

lossB = crossval(@(Xtrain,Ytrain,Xtest,Ytest) ...
    logLoss(glmfit(Xtrain, Ytrain, 'binomial'), Xtest, Ytest), ...
    sigmaMPoisson(:), confReports_Poiss(:), 'partition', cv);

fprintf("----------------------------------------\n " );
fprintf("Ground truth: Poisson decoder\n" );
fprintf("----------------------------------------\n" );
fprintf("Poisson decoder: %.4f \n", -mean(lossA) );
fprintf("Modulated Poisson decoder: %.4f \n", -mean(lossB) );

cvLL_GT_Poiss = [-lossA -lossB]; % Poiss, MPoiss

warning('on','all')

% how should i plot the results

% Modulated poisson - if no difference then can decoding happen?
%%

% Prop confidence by uncertainty level
confReports = {};

trlIdx = 1:ntrials; trlIdx = trlIdx(trlFltIdx);
trlCnt = sum(trlFltIdx);

for i = 1:trlCnt %1:ntrials
    key1 = sprintf('c_%g', trlMatrix(i, 1));
    key2 = sprintf('s_%g', trlMatrix(i, 2));
    key1 = matlab.lang.makeValidName(key1);
    key2 = matlab.lang.makeValidName(key2);

    confReports.MPoiss.(key1).(key2) = [];
    confReports.Poiss.(key1).(key2)  = [];
end

for i=1:trlCnt %1:ntrials
    key1 = sprintf('c_%g', trlMatrix(i, 1));
    key2 = sprintf('s_%g', trlMatrix(i, 2));
    key1 = matlab.lang.makeValidName(key1);
    key2 = matlab.lang.makeValidName(key2);
    
    confReports.MPoiss.(key1).(key2) = [confReports.MPoiss.(key1).(key2) confReports_MPoiss(i)];
    confReports.Poiss.(key1).(key2)  = [confReports.Poiss.(key1).(key2) confReports_Poiss(i)];
end

% Plotting business
pLC_MPoiss = [];
pLC_Poiss  = [];

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        pLC_MPoiss = [pLC_MPoiss 1 - mean( confReports.MPoiss.(key1).(key2) )];
        pLC_Poiss  = [pLC_Poiss 1 - mean( confReports.Poiss.(key1).(key2) )];
    end
end

figure
subplot(3, 3, 1)
plot(pLC_MPoiss, DisplayName="Modulated Poisson")
hold on
plot(pLC_Poiss, DisplayName="Poisson")
hold off
xlabel("Uncertainty conditions")
ylabel("P(LC)")
legend


% 
binEdges = linspace(min([sigmaPoisson sigmaMPoisson]), ...
    max([sigmaPoisson sigmaMPoisson]), 50);
binCenters = ( binEdges(1:end-1) + binEdges(2:end) ) / 2;

[binIdxMP, ~] = discretize(sigmaMPoisson, binEdges);
[binIdxP, ~]  = discretize(sigmaPoisson, binEdges);

T = table(confReports_MPoiss(:), confReports_Poiss(:), binIdxMP(:), binIdxP(:), ...
    'VariableNames', {'Conf_MP','Conf_P','BinID_MP','BinID_P'});

% Ground truth behavior: Modulated Poisson
confFn_MP_T = groupsummary(T, {'BinID_MP'}, {'mean'}, 'Conf_MP' ); % GT is Conf_MP always
confFn_P_T  = groupsummary(T, {'BinID_P'}, {'mean'}, 'Conf_MP' ); % GT is Conf_MP always

confFn_MP = zeros(1, numel(binCenters)) + nan;
confFn_P  = zeros(1, numel(binCenters)) + nan;

confFn_MP(confFn_MP_T.BinID_MP) = confFn_MP_T.mean_Conf_MP;
confFn_P(confFn_P_T.BinID_P(1:end-1)) = confFn_P_T.mean_Conf_MP(1:end-1);

confFn_MP = confFn_MP(~isnan(confFn_MP));
confFn_P  = confFn_P(~isnan(confFn_P));

subplot(3, 3, 2)
plot(binCenters(~isnan(confFn_MP)), 1 - confFn_MP, DisplayName="Modulated Poisson")
hold on
plot(binCenters(~isnan(confFn_P)), 1 - confFn_P, DisplayName="Poisson")
hold off
xlabel("\sigma_θ")
ylabel("Prop (LC)")
title("GT decoder: Modulated Poisson")
legend

% groupsummary(T, {'BinID_MP'}, {'mean'}, 'Conf_MP' )
% Ground truth behavior: Poisson 
confFn_MP_T = groupsummary(T, {'BinID_MP'}, {'mean'}, 'Conf_P' ); % GT is Conf_MP always
confFn_P_T  = groupsummary(T, {'BinID_P'}, {'mean'}, 'Conf_P' ); % GT is Conf_MP always

confFn_MP = zeros(1, numel(binCenters)) + nan;
confFn_P  = zeros(1, numel(binCenters)) + nan;

confFn_MP(confFn_MP_T.BinID_MP) = confFn_MP_T.mean_Conf_P;
confFn_P(confFn_P_T.BinID_P(1:end-1)) = confFn_P_T.mean_Conf_P(1:end-1);

confFn_MP = confFn_MP(~isnan(confFn_MP));
confFn_P  = confFn_P(~isnan(confFn_P));

subplot(3, 3, 3)
plot(binCenters(~isnan(confFn_MP)), 1 - confFn_MP, DisplayName="Modulated Poisson")
hold on
plot(binCenters(~isnan(confFn_P)), 1 - confFn_P, DisplayName="Poisson")
hold off
xlabel("\sigma_θ")
ylabel("Prop (LC)")
title("GT decoder: Poisson")
legend

colors = ['r', 'g', 'b', 'k', 'm'];  % define a color per point
subplot(3, 3, 4)
labels = {'Poisson dec', 'MPoiss dec'};
h = gobjects(2,1);
for i = 1:length(recoveryLL_GT_MPoiss)
    h(i)=stem(recoveryLL_GT_MPoiss(i), 1, Color=colors(i), LineWidth=1.5);
    hold on;
end
% stem(recoveryLL_GT_MPoiss, [1 1], Color=color)
xlabel("Log Likelihood")
legend(h, labels)
ylim([0 5])
xlim([-300 50])
title("GT: Modulated Poisson")

subplot(3, 3, 5)
h = gobjects(2,1);
for i = 1:length(recoveryLL_GT_Poiss)
    h(i) = stem(recoveryLL_GT_Poiss(i), 1, Color=colors(i), LineWidth=1.5);
    hold on;
end
xlabel("Log Likelihood")
ylim([0 5])
xlim([-300 50])
legend(h, labels)
title("GT: Poisson")

subplot(3, 3, 6)
% Ideally 
hold on
scatter([1 1 1 1 1], cvLL_GT_MPoiss(:,1), 70,'filled',  DisplayName='Poisson')
scatter([2 2 2 2 2], cvLL_GT_MPoiss(:,2), 70,'filled',  DisplayName='MPoisson')
hold off
hold on
xticks([1 2])
xticklabels(['Poisson', 'MPoiss'])
xlabel("Decoder")
ylabel("Log likelihood")
xlim([0 3])
legend
title(sprintf("Cross-validation\nGT: Modulated Poisson"))
% recoveryLL_GT_Poiss;
% cvLL_GT_Poiss

subplot(3, 3, 7)
% Ideally 
hold on
scatter([1 1 1 1 1], cvLL_GT_Poiss(:,1), 70,'filled',  DisplayName='Poisson')
scatter([2 2 2 2 2], cvLL_GT_Poiss(:,2),  70, 'filled', DisplayName='MPoisson')
hold off
hold on
xticks([1 2])
xticklabels(['Poisson', 'MPoiss'])
xlabel("Decoder")
ylabel("Log likelihood")
legend
xlim([0 3])
title(sprintf("Cross-validation\nGT: Poisson"))
% recoveryLL_GT_Poiss;
% cvLL_GT_Poiss

%%
figure

colors = {};
itr = 1;
for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        if itr == 1
            colors.(key1).(key2) = [1 0.6 0.6]; % light red 
        elseif itr == 2
            colors.(key1).(key2) = [0.6 0.6 1]; % light blue
        elseif itr == 3
            colors.(key1).(key2) = [1 0 0]; % red 
        elseif itr == 4
            colors.(key1).(key2) = [0 0 1]; % blue
        end

        itr = itr + 1;
    end
end

% tuning functions
subplot(2, 4, 1)
hold on

nrnIdx = 69; %randi(nNeurons); % Randomly select one neuron index, 69 is a nice index

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        nrnTuningFn = tuningFnsByExpConditions.(key1).(key2)(nrnIdx, :);
        
        plot(rad2deg(tuningFnOriSpace), nrnTuningFn, ...
            DisplayName=sprintf("(%.2f, %d)", contrasts(cIdx), spreads(sIdx)), ...
            LineWidth=2, Color=colors.(key1).(key2))
        xlabel("Orientation (deg)")
        ylabel("Response (IPS)")
        legend
        title("1 unit")
    end
end

hold off

% Mean vs variance plot
subplot(2, 4, 2)
hold on

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        m = squeeze(meanSpkCnt(nrnIdx, cIdx, sIdx, :));
        v = squeeze(varSpkCnt(nrnIdx, cIdx, sIdx, :));

        % Define custom model (edit as needed)
        x = m; y = v;
        ft = fittype('x + sigma_g^2*(x)^2', ...
                     'independent','x','coefficients',{'sigma_g'});
        
        opts = fitoptions(ft);
        opts.StartPoint = [rand];     % initial guess
        opts.Lower = [0];          % e.g., sigma_g >= 0
        % opts.Upper = [10];       % optional upper bound
        
        [curve, ~] = fit(x, y, ft, opts);
        coeffs = coeffvalues(curve);
        
        xFit = linspace(0, max(x), 100);
        yFit = xFit + coeffs^2*(xFit).^2;

        scatter(m, v, "filled", "MarkerFaceAlpha", 0.3, ...
            "MarkerFaceColor", colors.(key1).(key2), ...
            "MarkerEdgeColor", colors.(key1).(key2), ...
            DisplayName=sprintf("(%.2f, %d)", contrasts(cIdx), spreads(sIdx)), ...
            Color=colors.(key1).(key2))
        plot([min(m) max(m)], [min(m) max(m)], 'k--', HandleVisibility='off') 
        plot(xFit, yFit, HandleVisibility="off", LineWidth=2, Color=colors.(key1).(key2))
        
        xlabel("Mean Spk Cnt")
        ylabel("Var Spk Cnt")
        legend
        title("1 unit")
    end
end

hold off

% Likelihood - Modulated Poisson at different condition
subplot(2, 4, 3)
hold on

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        trlPdf = pdfOriByCondition.MPoisson.(key1).(key2).o_1_91986;
        ori = rad2deg(1.91986);
        
        plot(rad2deg(tuningFnOriSpace), trlPdf, LineWidth=2, Color=colors.(key1).(key2))
        
        xlabel("Orientation (deg)")
        ylabel("P(ori)")
        
        xStart = ori - 15;
        xEnd   = ori + 15;
        xlim([xStart, xEnd])
        title("Modulated Poisson (1 trial)")

    end
end

hold off

% Likelihood - Poisson at different condition
subplot(2, 4, 4)
hold on

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        trlPdf = pdfOriByCondition.Poisson.(key1).(key2).o_1_91986;
        ori = rad2deg(1.91986);
        
        plot(rad2deg(tuningFnOriSpace), trlPdf, LineWidth=2, Color=colors.(key1).(key2))
        
        xlabel("Orientation (deg)")
        ylabel("P(ori)")
        
        xStart = ori - 15;
        xEnd   = ori + 15;
        xlim([xStart, xEnd])
        title("Poisson (1 trial)")

    end
end

hold off

% Estimation error : MPoisson
subplot(2, 4, 5)
hold on

itr = 1;
for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        err = estimationErr.MPoisson.(key1).(key2);

        xPts = itr + 0.4*(rand(1, numel(err)) - 0.5);
        
        scatter(xPts, err, 'filled', "MarkerFaceAlpha", 0.3, ...
            "MarkerFaceColor", colors.(key1).(key2), ...
            "MarkerEdgeColor", colors.(key1).(key2));    
        plot([itr - 0.35, itr + 0.35], [mean(err), mean(err)], ...
            LineStyle="-", LineWidth=2, Color='black')
        
        ylabel("Estimation err (deg)")
        ylim([-12 12])
        title("Modulated Poisson")

        itr = itr + 1;
    end
end
hold off

% Estimation error : Poisson
subplot(2, 4, 6)
hold on

itr = 1;
for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        err = estimationErr.Poisson.(key1).(key2);
        
        xPts = itr + 0.4*(rand(1, numel(err)) - 0.5);
        
        scatter(xPts, err, 'filled', "MarkerFaceAlpha", 0.3, ...
            "MarkerFaceColor", colors.(key1).(key2), ...
            "MarkerEdgeColor", colors.(key1).(key2));    
        plot([itr - 0.35, itr + 0.35], [mean(err), mean(err)], ...
            LineStyle="-", LineWidth=2, Color='black')
        
        ylabel("Estimation err (deg)")
        ylim([-12 12])
        title("Poisson")

        itr = itr + 1;
    end
end
hold off

subplot(2, 4, 7)
hold on

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        err = estimationErr.MPoisson.(key1).(key2);
        sigmas = rad2deg( sigma_theta.MPoisson.(key1).(key2) );
        
        scatter(sigmas, err, 'filled', "MarkerFaceAlpha", 0.3, ...
            "MarkerFaceColor", colors.(key1).(key2), ...
            "MarkerEdgeColor", colors.(key1).(key2));  

        xlabel("\sigma_θ (deg)")
        ylabel("Estimation err (deg)")
        title("Modulated Poisson")
    end
end

hold off

subplot(2, 4, 8)
hold on

for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        err = estimationErr.Poisson.(key1).(key2);
        sigmas = rad2deg( sigma_theta.Poisson.(key1).(key2) );
        
        scatter(sigmas, err, 'filled', "MarkerFaceAlpha", 0.3, ...
            "MarkerFaceColor", colors.(key1).(key2), ...
            "MarkerEdgeColor", colors.(key1).(key2));  

        xlabel("\sigma_θ (deg)")
        ylabel("Estimation err (deg)")
        title("Poisson")
    end
end

hold off

%% Figure2
figure

colorVector = zeros(ntrials, 3);
sigmaPoisson = zeros(ntrials, 1);
sigmaMPoisson = zeros(ntrials, 1);

for i = 1:ntrials
    key1 = sprintf('c_%g', trlMatrix(i, 1));
    key2 = sprintf('s_%g', trlMatrix(i, 2));
    key1 = matlab.lang.makeValidName(key1);
    key2 = matlab.lang.makeValidName(key2);

    colorVector(i, :) = colors.(key1).(key2);

    sigmaPoisson(i)  = tuningFnData.data.trialData.Poisson.metrics{i}.sigma;
    sigmaMPoisson(i) = tuningFnData.data.trialData.MPoisson.metrics{i}.sigma;
end

subplot(2, 3, 1)
scatter(decodedThetasPossDec, decodedThetasMPossDec, 36, colorVector, ...
    'filled', "MarkerFaceAlpha", 0.7);  
xlabel("Decoded theta (Poisson)")
ylabel("Decoded theta (Modulated Poisson)")

subplot(2, 3, 2)
scatter(decodedContrastsPossDec, decodedContrastsMPossDec, 36, colorVector, ...
    'filled', "MarkerFaceAlpha", 0.3);  
xlabel("Decoded contrast (Poisson)")
ylabel("Decoded contrast (Modulated Poisson)")

subplot(2, 3, 3)
scatter(decodedSpreadsPossDec, decodedSpreadsMPossDec, 36, colorVector, ...
    'filled', "MarkerFaceAlpha", 0.3);  
xlabel("Decoded spread (Poisson)")
ylabel("Decoded spread (Modulated Poisson)")

subplot(2, 3, 4)
scatter(sigmaPoisson, sigmaMPoisson, 36, colorVector, ...
    'filled', "MarkerFaceAlpha", 0.3);  
hold on
plot([0 max(sigmaPoisson)], [0 max(sigmaPoisson)], 'k--', HandleVisibility='off')
hold off
xlabel("\sigma_θ (Poisson)")
ylabel("\sigma_θ (Modulated Poisson)")

subplot(2, 3, 5)
hold on

itr = 1;
for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);

        err = decodedThetaDiff.(key1).(key2);
        
        xPts = itr + 0.4*(rand(1, numel(err)) - 0.5);
        
        scatter(xPts, err, 'filled', "MarkerFaceAlpha", 0.3, ...
            "MarkerFaceColor", colors.(key1).(key2), ...
            "MarkerEdgeColor", colors.(key1).(key2));    
        
        ylabel("θ (MPoisson - Poisson)")
        
        itr = itr + 1;
    end
end
hold off

subplot(2, 3, 6)
hold on

itr = 1;
for cIdx = 1:numel(contrasts)
    for sIdx = 1:numel(spreads)
        key1 = sprintf('c_%g', contrasts(cIdx));
        key2 = sprintf('s_%g', spreads(sIdx));
        key1 = matlab.lang.makeValidName(key1);
        key2 = matlab.lang.makeValidName(key2);
        
        err = decodedSigmaDiff.(key1).(key2);
        
        xPts = itr + 0.4*(rand(1, numel(err)) - 0.5);
        
        scatter(xPts, err, 'filled', "MarkerFaceAlpha", 0.3, ...
            "MarkerFaceColor", colors.(key1).(key2), ...
            "MarkerEdgeColor", colors.(key1).(key2));    
        
        ylabel("\sigma_θ (MPoisson - Poisson)")
        
        itr = itr + 1;
    end
end
hold off

%%
% function varGain = getVarGain(spread, contrast)
%     varGain = 0.0001 * spread ./ contrast;
% end

% getVarGain(30, 0.05)

function varGain = getVarGain(spread, contrast)
%     sigmaG = (0.15/27)*spread + (0.25 - (0.15/27)*30 ) + ...
%         ( - (0.2/0.04)*contrast + (0.3 + (0.2/0.04)*0.01) + 0.4);
    sigmaG = (0.15/27)*spread + (0.25 - (0.15/27)*30 ) + ...
        ( - (0.2/0.04)*contrast + (0.3 + (0.2/0.04)*0.01) + 0.4);
    varGain = sigmaG.^2;
end

function L = logLoss(b, Xtest, Ytest)
    p = glmval(b, Xtest, 'logit');
    eps = 1e-10;
    L = -sum(Ytest.*log(p+eps) + (1-Ytest).*log(1-p+eps));
end