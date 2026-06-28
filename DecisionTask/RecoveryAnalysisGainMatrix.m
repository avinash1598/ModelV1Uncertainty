clear all
%close all
clc

rng('shuffle');

% TuningFn on remote computer

addpath('C:\Users\avinash1598\Desktop\ModelV1Uncertainty\DecisionTask\Scripts\')
%addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/ModelV1Uncertainty/EstimationTask/Scripts/')
% addpath('/Volumes/Desktop/UT Austin/Goris lab/ModelV1Uncertainty/DecisionTask')

%tuningFnData = load('tuningFnDataGainMatrix.mat'); % Decode using gains of individual neurons
%tuningFnData = load('tuningFnDataDecodeAvgGain.mat'); % Decode using average gains of all the neurons
tuningFnData = load('tuningFnDataCovSigmaG_v2.mat'); % Decode using average gains of all the neurons

% How is average decoder compared to the original decoder (individual gains
% for each neuron)?
% Which decoder (avg or full) better explain the actual behavioral results from monkey? 

% ----------------------------------
% Params
% ----------------------------------
nNeurons     = 200;   % Count of neurons
stimDuration = 2;     % Stimulus duration set to 1 seconds
timeStep     = 0.001; % 0.001s (1ms) - Step size of time bins used for binning stimulus duration 

% Stimulus parameters (angles in radians)
warning("don't change these values")
contrasts               = [0.01 0.05]; 
spreads                 = deg2rad([3 30]); 
uniqStimOris            = linspace(90-10, 90+10, 21); %80:0.5:100;
uniqStimOris            = deg2rad(uniqStimOris');
stimParam.numStim       = numel(uniqStimOris);                                % Number of unique stimuli
stimParam.countPerStim  = 60; % 100 
ntrials                 = stimParam.numStim * stimParam.countPerStim * numel(contrasts) * numel(spreads);  % Total number of trials

trlMatrix               = tuningFnData.data.trialMatrix;
%%
timeBins                  = 0:timeStep:stimDuration; 
stimRespProfile           = 1 + zeros(1, numel(timeBins));
%gainVector                = gainVector'; %nNeurons x nTrials 1 + zeros(nNeurons, ntrials); % constant gain - NO gain modulation

tuningParams              = tuningFnData.data.tuningParams;
neuronsPrefOrientation    = tuningFnData.data.neuronsPrefOrientation;

tuningFnOriSpace  = linspace(0, pi, 361);
tuningFnContrasts = linspace(1e-4, 0.15, 49);
tuningFnSpreads   = deg2rad(linspace(1, 90, 50));

% Update tuning Fns
tuningFns    = tuningFnData.data.tuningFns;
nrnVarGains  = tuningFnData.data.nrnVarGains;

%% LL comparison
warning('off')

% Is decoded uncertainty important for predicting confidence?
% Signal strength alone (seems like this alone is not sufficient for explaining confidence)
% Decoded uncertainty alone
% Combined confidence variable

% LL comparison 1: decision variable: Poisson vs M poisson
% LL comparison 2: sigma : Poisson vs M poisson
% LL comparison 3: confidence variable : Poisson vs M poisson

sigmaPoisson = [];
sigmaMPoisson = [];
for i=1:ntrials
    sigmaPoisson = [sigmaPoisson tuningFnData.data.trialData.Poisson.metrics{i}.sigma];
    sigmaMPoisson = [sigmaMPoisson tuningFnData.data.trialData.MPoisson.metrics{i}.sigma];
end

decisionVarPoiss  = tuningFnData.data.decodedThetasPossDec - deg2rad(90);
decisionVarMPoiss = tuningFnData.data.decodedThetasMPossDec - deg2rad(90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ground truth: Poisson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get decision and confidence from poisson decoder
confVarPoissDec    = tuningFnData.data.confVarPoissDec;
p = 50; Cc_Poiss   = prctile(confVarPoissDec, p);

confReportsPoissDec  = confVarPoissDec > Cc_Poiss;
decisionPoissDec     = tuningFnData.data.decodedThetasPossDec >= deg2rad(90);

% Get confidence variable from poisson and modulated poisson decoder
confVarPoissDec    = tuningFnData.data.confVarPoissDec;
confVarMPoissDec   = tuningFnData.data.confVarMPoissDec;

figure

for i = 1:4

    % GT: Poisson
    if i == 1
        % Logistic regression on decision variable
        outputLabels  = decisionPoissDec(:); % GT: Poisson
        predictors    = decisionVarPoiss; % decision var

        x_label = 'Decision var';
        y_label = 'Decision';

        title_txt = 'Choice recovery';
    elseif i == 2
        % Logistic regression on confidence variable
        outputLabels  = confReportsPoissDec(:); % GT: Poisson
        predictors    = decisionVarPoiss; % decision var

        x_label = 'Decision var';
        y_label = 'Confidence';

        title_txt = 'Confidence recovery';
    elseif i == 3
        % Logistic regression on confidence variable
        outputLabels  = confReportsPoissDec(:); % GT: Poisson
        predictors    = sigmaPoisson; % decoded sigma

        x_label = 'Sigma';
        y_label = 'Confidence';

        title_txt = 'Confidence recovery';
    else
        outputLabels  = confReportsPoissDec(:); % GT: Poisson
        predictors    = confVarPoissDec; % confidence variable

        x_label = 'Confidence variable';
        y_label = 'Confidence';

        title_txt = 'Confidence recovery';
    end

    [ll_Poiss, xFit, pFit] = computeLikelihood(outputLabels, predictors);
    losses_Poiss           = crossValidate(outputLabels, predictors);
    
    subplot(3, 4, i)
    scatter(predictors, outputLabels, DisplayName='Data'); hold on
    plot(xFit, pFit, DisplayName='Logistic fit'); hold off
    xlabel(x_label)
    ylabel(y_label)
    legend
    title('Poisson')

    % M Poisson
    if i == 1
        % Logistic regression on decision variable
        outputLabels  = decisionPoissDec(:); % GT: Poisson
        predictors    = decisionVarMPoiss; % decision var

        x_label = 'Decision var';
        y_label = 'Decision';
    elseif i == 2
        % Logistic regression on confidence variable
        outputLabels  = confReportsPoissDec(:); % GT: Poisson
        predictors    = decisionVarMPoiss; % decision var

        x_label = 'Decision var';
        y_label = 'Confidence';
    elseif i == 3
        % Logistic regression on confidence variable
        outputLabels  = confReportsPoissDec(:); % GT: Poisson
        predictors    = sigmaMPoisson; % decoded sigma

        x_label = 'Sigma';
        y_label = 'Confidence';
    else
        outputLabels  = confReportsPoissDec(:); % GT: Poisson
        predictors    = confVarMPoissDec; % confidence variable

        x_label = 'Confidence variable';
        y_label = 'Confidence';
    end
    
    [ll_MPoiss, xFit, pFit] = computeLikelihood(outputLabels, predictors);
    losses_MPoiss           = crossValidate(outputLabels, predictors);
    
    subplot(3, 4, 4 + i)
    scatter(predictors, outputLabels, DisplayName='Data'); hold on
    plot(xFit, pFit, DisplayName='Logistic fit'); hold off
    xlabel(x_label)
    ylabel(y_label)
    legend
    title('M Poisson')

    colors = ['r', 'g', 'b', 'k', 'm'];  % define a color per point
    labels = {'Poisson dec', 'MPoiss dec'};

    subplot(3, 4, 8 + i)
    h = gobjects(2,1);
    h(1) = stem(ll_Poiss, 1, Color=colors(1), LineWidth=1.5); hold on
    h(2) = stem(ll_MPoiss, 1, Color=colors(2), LineWidth=1.5); hold off
    hold on;
    xlabel("Log Likelihood")
    ylim([0 5])
    
    minV = min([ll_Poiss ll_MPoiss]);
    rangeV = max([ll_Poiss ll_MPoiss]) - min([ll_Poiss ll_MPoiss]);
    delta = (2000 - rangeV);
    
    xlim([minV - 10 - delta / 2, minV + rangeV + delta / 2])
    %xlim(minV + [-900 10])
    legend(h, labels)
    title("GT: Poisson (" + title_txt + ")")
    
    fprintf(title_txt + "\n")
    fprintf("Log Likelihood comparison & CV\n\n")
    fprintf("----------------------------------------\n " );
    fprintf("Ground truth: Poisson decoder\n" );
    fprintf("----------------------------------------\n" );
    fprintf("Poisson decoder: %.4f \n", ll_Poiss );
    fprintf("Modulated Poisson decoder: %.4f \n", ll_MPoiss );
    fprintf("CV: Poisson decoder: %.4f \n", mean(losses_Poiss) );
    fprintf("CV: Modulated Poisson decoder: %.4f \n", mean(losses_MPoiss) );
    fprintf("\n")
end


% Compare decision variace strength predictability as well

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ground truth: Modulated Poisson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get decision and confidence from poisson decoder
confVarMPoissDec    = tuningFnData.data.confVarMPoissDec;
p = 50; Cc_MPoiss   = prctile(confVarMPoissDec, p);

confReportsMPoissDec  = confVarMPoissDec > Cc_MPoiss;
decisionMPoissDec     = tuningFnData.data.decodedThetasMPossDec >= deg2rad(90);

% Get confidence variable from poisson and modulated poisson decoder
confVarPoissDec    = tuningFnData.data.confVarPoissDec;
confVarMPoissDec   = tuningFnData.data.confVarMPoissDec;

figure

for i = 1:4

    % Poisson
    if i == 1
        % Logistic regression on decision variable
        outputLabels  = decisionMPoissDec(:); % GT: M Poisson
        predictors    = decisionVarPoiss; % decision var

        x_label = 'Decision var';
        y_label = 'Decision';

        title_txt = 'Choice recovery';
    elseif i == 2
        % Logistic regression on confidence variable
        outputLabels  = confReportsMPoissDec(:); % GT: M Poisson
        predictors    = decisionVarPoiss; % decision var

        x_label = 'Decision var';
        y_label = 'Confidence';

        title_txt = 'Confidence recovery';
    elseif i == 3
        % Logistic regression on confidence variable
        outputLabels  = confReportsMPoissDec(:); % GT: M Poisson
        predictors    = sigmaPoisson; % decoded sigma

        x_label = 'Sigma';
        y_label = 'Confidence';

        title_txt = 'Confidence recovery';
    else
        outputLabels  = confReportsMPoissDec(:); % GT: M Poisson
        predictors    = confVarPoissDec; % confidence variable

        x_label = 'Confidence variable';
        y_label = 'Confidence';

        title_txt = 'Confidence recovery';
    end
    
    [ll_Poiss, xFit, pFit] = computeLikelihood(outputLabels, predictors);
    losses_Poiss           = crossValidate(outputLabels, predictors);
    
    subplot(3, 4, i)
    scatter(predictors, outputLabels, DisplayName='Data'); hold on
    plot(xFit, pFit, DisplayName='Logistic fit'); hold off
    xlabel(x_label)
    ylabel(y_label)
    legend
    title('Poisson')
    
    % M Poisson
    if i == 1
        % Logistic regression on decision variable
        outputLabels  = decisionMPoissDec(:); % GT: M Poisson
        predictors    = decisionVarMPoiss; % decision var
        
        x_label = 'Decision var';
        y_label = 'Decision';
    elseif i == 2
        % Logistic regression on confidence variable
        outputLabels  = confReportsMPoissDec(:); % GT: M Poisson
        predictors    = decisionVarMPoiss; % decision var

        x_label = 'Decision var';
        y_label = 'Confidence';
    elseif i == 3
        % Logistic regression on confidence variable
        outputLabels  = confReportsMPoissDec(:); % GT: M Poisson
        predictors    = sigmaMPoisson; % decoded sigma

        x_label = 'Sigma';
        y_label = 'Confidence';
    else
        outputLabels  = confReportsMPoissDec(:); % GT: M Poisson
        predictors    = confVarMPoissDec; % confidence variable

        x_label = 'Confidence variable';
        y_label = 'Confidence';
    end
    
    [ll_MPoiss, xFit, pFit] = computeLikelihood(outputLabels, predictors);
    losses_MPoiss           = crossValidate(outputLabels, predictors);
    
    subplot(3, 4, 4 + i)
    scatter(predictors, outputLabels, DisplayName='Data'); hold on
    plot(xFit, pFit, DisplayName='Logistic fit'); hold off
    xlabel(x_label)
    ylabel(y_label)
    legend
    title('M Poisson')

    colors = ['r', 'g', 'b', 'k', 'm'];  % define a color per point
    labels = {'Poisson dec', 'MPoiss dec'};

    subplot(3, 4, 8 + i)
    h = gobjects(2,1);
    h(1) = stem(ll_Poiss, 1, Color=colors(1), LineWidth=1.5); hold on
    h(2) = stem(ll_MPoiss, 1, Color=colors(2), LineWidth=1.5); hold off
    hold on;
    xlabel("Log Likelihood")
    ylim([0 5])
    
    minV = min([ll_Poiss ll_MPoiss]);
    rangeV = max([ll_Poiss ll_MPoiss]) - min([ll_Poiss ll_MPoiss]);
    delta = (2000 - rangeV);
    
    xlim([minV - 10 - delta / 2, minV + rangeV + delta / 2])
    %xlim(minV + [-900 10])
    legend(h, labels)
    title("GT: M Poisson (" + title_txt + ")")
    
    fprintf(title_txt + "\n")
    fprintf("Log Likelihood comparison\n\n")
    fprintf("----------------------------------------\n " );
    fprintf("Ground truth: Modulated Poisson decoder\n" );
    fprintf("----------------------------------------\n" );
    fprintf("Poisson decoder: %.4f \n", ll_Poiss );
    fprintf("Modulated Poisson decoder: %.4f \n", ll_MPoiss );
    fprintf("CV: Poisson decoder: %.4f \n", mean(losses_Poiss) );
    fprintf("CV: Modulated Poisson decoder: %.4f \n", mean(losses_MPoiss) );
    fprintf("\n")
end


warning('on')

% TODO: plot distributions of LLs
% confidence: decoded uncertainty which decoder 
%% Psychometric functions
% close all

decisionVarPoiss   = rad2deg(tuningFnData.data.decodedThetasPossDec);
decisionPoissDec   = tuningFnData.data.decodedThetasPossDec >= deg2rad(90);
confVarPoissDec    = tuningFnData.data.confVarPoissDec;
decisionVarMPoiss  = rad2deg(tuningFnData.data.decodedThetasMPossDec);
decisionMPoissDec  = tuningFnData.data.decodedThetasMPossDec >= deg2rad(90);
confVarMPoissDec   = tuningFnData.data.confVarMPoissDec;

p = 50;
Cc_Poiss  = prctile(confVarPoissDec, p);
Cc_MPoiss = prctile(confVarMPoissDec, p);

confPoiss  = confVarPoissDec > Cc_Poiss; % HC: 1, LC: 0
confMPoiss = confVarMPoissDec > Cc_MPoiss; % HC: 1, LC: 0

% decodedSigmaPoiss  = tuningFnData.data.decodedSpreadsPossDec;
% decodedSigmaMPoiss = tuningFnData.data.decodedSpreadsMPossDec;
% decodedThetaPoiss  = tuningFnData.data.decodedSpreadsPossDec;
% decodedThetaMPoiss = tuningFnData.data.decodedSpreadsMPossDec;

% % 4 conditions
% psycFn_HC_GTPoiss = zeros(4, stimParam.numStim);
% psycFn_LC_GTPoiss = zeros(4, stimParam.numStim);
% 
% psycFn_HC_GTMPoiss = zeros(4, stimParam.numStim);
% psycFn_LC_GTMPoiss = zeros(4, stimParam.numStim);

psycFn_HC_GTPoiss = zeros(1, stimParam.numStim);
psycFn_LC_GTPoiss = zeros(1, stimParam.numStim);
psycFn_HC_GTMPoiss = zeros(1, stimParam.numStim);
psycFn_LC_GTMPoiss = zeros(1, stimParam.numStim);

countFn_HC_GTPoiss = zeros(1, stimParam.numStim);
countFn_LC_GTPoiss = zeros(1, stimParam.numStim);
countFn_HC_GTMPoiss = zeros(1, stimParam.numStim);
countFn_LC_GTMPoiss = zeros(1, stimParam.numStim);

confFn_HC_GTPoiss  = zeros(1, stimParam.numStim);
confFn_LC_GTPoiss  = zeros(1, stimParam.numStim);
confFn_HC_GTMPoiss = zeros(1, stimParam.numStim);
confFn_LC_GTMPoiss = zeros(1, stimParam.numStim);

trlContrasts = trlMatrix(:, 1);
trlSpreads   = trlMatrix(:, 2);
trlOris      = trlMatrix(:, 3);

idx = 1;
% for cIdx = 1:numel(contrasts)
%     for sIdx = 1:numel(spreads)
%         key1 = sprintf('c_%g', contrasts(cIdx));
%         key2 = sprintf('s_%g', spreads(sIdx));
%         key1 = matlab.lang.makeValidName(key1);
%         key2 = matlab.lang.makeValidName(key2);

        for oriIdx = 1:numel(uniqStimOris)
            
%             fltIdx = ( trlContrasts == contrasts(cIdx) ) & ...
%                 ( trlSpreads == spreads(sIdx) ) & ...
%                 ( trlOris == uniqStimOris(oriIdx) );

            fltIdx = ( trlOris == uniqStimOris(oriIdx) );
            
            % Poisson
            fltDecisions_Poiss = decisionPoissDec(fltIdx);
            fltConf_Poiss      = confPoiss(fltIdx);
            
            psycFn_HC_GTPoiss(idx, oriIdx) = sum( (fltDecisions_Poiss == 1) & (fltConf_Poiss == 1)  ) ...  % & (fltConf_Poiss == 1) 
                / numel(fltDecisions_Poiss((fltConf_Poiss == 1))); % Prop CCW (or is it CW? doesn't matter though)
            psycFn_LC_GTPoiss(idx, oriIdx) = sum( (fltDecisions_Poiss == 1) & (fltConf_Poiss == 0) ) ... %  & (fltConf_Poiss == 0)
                / numel(fltDecisions_Poiss((fltConf_Poiss == 0))); % Prop CCW (or is it CW? doesn't matter though)
            
            % Count of HC and LC at each orientation
            %countFn_HC_GTPoiss(idx, oriIdx) = sum( (fltDecisions_Poiss == 1) & (fltConf_Poiss == 1)  );
            %countFn_LC_GTPoiss(idx, oriIdx) = sum( (fltDecisions_Poiss == 1) & (fltConf_Poiss == 0)  );
            countFn_HC_GTPoiss(idx, oriIdx) = sum( (fltConf_Poiss == 1)  );
            countFn_LC_GTPoiss(idx, oriIdx) = sum( (fltConf_Poiss == 0)  );
            
            % Proportion HC and LC at each stim
            confFn_HC_GTPoiss(idx, oriIdx) = numel(fltDecisions_Poiss((fltConf_Poiss == 1))) / numel(fltDecisions_Poiss);
            confFn_LC_GTPoiss(idx, oriIdx) = numel(fltDecisions_Poiss((fltConf_Poiss == 0))) / numel(fltDecisions_Poiss);
            
            % Modulated Poisson
            fltDecisions_MPoiss = decisionMPoissDec(fltIdx);
            fltConf_MPoiss      = confMPoiss(fltIdx);
            
            psycFn_HC_GTMPoiss(idx, oriIdx) = sum( (fltDecisions_MPoiss == 1) & (fltConf_MPoiss == 1)  ) ... % & (fltConf_MPoiss == 1)
                / numel(fltDecisions_MPoiss(fltConf_MPoiss == 1)); % Prop CCW (or is it CW? doesn't matter though)
            psycFn_LC_GTMPoiss(idx, oriIdx) = sum( (fltDecisions_MPoiss == 1) & (fltConf_MPoiss == 0) ) ... % & (fltConf_MPoiss == 0)
                / numel(fltDecisions_MPoiss(fltConf_MPoiss == 0)); % Prop CCW (or is it CW? doesn't matter though)
            
            %countFn_HC_GTMPoiss(idx, oriIdx) = sum( (fltDecisions_MPoiss == 1) & (fltConf_MPoiss == 1)  );
            %countFn_LC_GTMPoiss(idx, oriIdx) = sum( (fltDecisions_MPoiss == 1) & (fltConf_MPoiss == 0)  );
            countFn_HC_GTMPoiss(idx, oriIdx) = sum( (fltConf_MPoiss == 1)  );
            countFn_LC_GTMPoiss(idx, oriIdx) = sum( (fltConf_MPoiss == 0)  );
            
            % Proportion HC and LC at each stim
            confFn_HC_GTMPoiss(idx, oriIdx) = numel(fltDecisions_MPoiss((fltConf_MPoiss == 1))) / numel(fltDecisions_MPoiss);
            confFn_LC_GTMPoiss(idx, oriIdx) = numel(fltDecisions_MPoiss((fltConf_MPoiss == 0))) / numel(fltDecisions_MPoiss);
            
        end
        
%         idx = idx + 1;
%     end
% end

figure

colors_highConf = [0.12 0.50 0.30];   % deep green
colors_lowConf  = [0.85 0.75 0.20];   % olive-yellow

% idx = 1;
% for cIdx = 1:numel(contrasts)
%     for sIdx = 1:numel(spreads)
        % Poisson
        subplot(2, 4, idx)
        plot(rad2deg(uniqStimOris), psycFn_HC_GTPoiss, DisplayName="HC", LineWidth=2, Color=colors_highConf) % (idx, :)
        hold on
        scatter(rad2deg(uniqStimOris), psycFn_HC_GTPoiss, 1 + countFn_HC_GTPoiss, ...
            'filled', 'MarkerFaceColor', colors_highConf, 'MarkerEdgeColor', 'w', HandleVisibility='off')

        plot(rad2deg(uniqStimOris), psycFn_LC_GTPoiss, DisplayName="LC", LineWidth=2, Color=colors_lowConf)
        scatter(rad2deg(uniqStimOris), psycFn_LC_GTPoiss, 1 + countFn_LC_GTPoiss, ...
            'filled', 'MarkerFaceColor', colors_lowConf, 'MarkerEdgeColor', 'w', HandleVisibility='off')
        xline(90, LineStyle="--", HandleVisibility='off')
        hold off
        xlabel("Orientation (deg)")
        ylabel("Proportion CCW")
        legend
        %title(sprintf("Poisson C: %.2f, S: %d", contrasts(cIdx), spreads(sIdx)))
        title(sprintf("Poisson"))
        
        set(gca, 'FontSize', 16, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'off')

        % Modulated Poisson
        subplot(2, 4, 4+idx)
        plot(rad2deg(uniqStimOris), psycFn_HC_GTMPoiss, DisplayName="HC", LineWidth=2, Color=colors_highConf)
        hold on
        scatter(rad2deg(uniqStimOris), psycFn_HC_GTMPoiss, 1 + countFn_HC_GTMPoiss, ...
            'filled', 'MarkerFaceColor', colors_highConf, 'MarkerEdgeColor', 'w', HandleVisibility='off')
        
        plot(rad2deg(uniqStimOris), psycFn_LC_GTMPoiss, DisplayName="LC", LineWidth=2, Color=colors_lowConf)
        scatter(rad2deg(uniqStimOris), psycFn_LC_GTMPoiss, 1 + countFn_LC_GTMPoiss, ...
            'filled', 'MarkerFaceColor', colors_lowConf, 'MarkerEdgeColor', 'w', HandleVisibility='off')
        xline(90, LineStyle="--", HandleVisibility='off')
        hold off
        xlabel("Orientation (deg)")
        ylabel("Proportion CCW")
        legend
        %title(sprintf("Modulated Poisson C: %.2f, S: %d", contrasts(cIdx), spreads(sIdx)))
        title(sprintf("Modulated Poisson"))
        
        set(gca, 'FontSize', 16, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'off')

        idx = idx + 1;
% 
%     end
% end

subplot(2, 4, 2)
plot(rad2deg(uniqStimOris), confFn_HC_GTPoiss)
xlabel("Ori (deg)")
ylabel("Prop HC")
title("Poisson")
ylim([0 1])

subplot(2, 4, 6)
plot(rad2deg(uniqStimOris), confFn_HC_GTMPoiss)
xlabel("Ori (deg)")
ylabel("Prop HC")
title("M Poisson")
ylim([0 1])

subplot(2, 4, 3)
scatter(decisionVarPoiss, confVarPoissDec, 'filled')
ylabel("Conf variable")
xlabel("Decision variable")
title("Poisson")

subplot(2, 4, 7)
scatter(decisionVarMPoiss, confVarMPoissDec, 'filled')
ylabel("Conf variable")
xlabel("Decision variable")
title("M Poisson")


totalTrlCnt = numel(decisionPoissDec);

figure
subplot(2, 3, 1)
scatter(rad2deg(tuningFnData.data.decodedThetasPossDec), rad2deg(tuningFnData.data.decodedThetasMPossDec))
xlabel("Theta (Poss decoder)")
ylabel("Theta (M Poss decoder)")

subplot(2, 3, 2)
scatter(decisionPoissDec + 0.1*rand(totalTrlCnt, 1), ...
    decisionMPoissDec + 0.1*rand(totalTrlCnt, 1), 'filled', 'MarkerFaceAlpha', 0.3)
xlabel("Decision (Poss decoder)")
ylabel("Decision (M Poss decoder)")

subplot(2, 3, 3)
scatter(confPoiss + 0.1*rand(totalTrlCnt, 1), ...
    confMPoiss + 0.1*rand(totalTrlCnt, 1), 'filled', 'MarkerFaceAlpha', 0.3)
xlabel("Confidence (Poss decoder)")
ylabel("Confidence (M Poss decoder)")

subplot(2, 3, 4)
scatter(confVarPoissDec, ...
    confVarMPoissDec, 'filled', 'MarkerFaceAlpha', 0.3)
hold on
plot([0, 35], [0, 35], LineStyle="--")
hold off
xlabel("Confidence Var (Poss)")
ylabel("Confidence Var (M Poss)")

subplot(2, 3, 5)
x = rad2deg(tuningFnData.data.decodedThetasPossDec - deg2rad(90));
y = rad2deg(tuningFnData.data.decodedThetasMPossDec - deg2rad(90));

scatter(x, ...
    y, 'filled', 'MarkerFaceAlpha', 0.3)
hold on
plot([0, 25], [0, 25], LineStyle="--")
hold off
xlabel("Decision Var (Poss)")
ylabel("Decision Var (M Poss)")

subplot(2, 3, 6)
scatter(sigmaPoisson, sigmaMPoisson, 'filled')
hold on
plot([0, max(sigmaMPoisson)], [0, max(sigmaMPoisson)], LineStyle="--")
hold off
xlabel("Sigma Poisson")
ylabel("Sigma M Poisson")


%% Functions
function varGain = getVarGain(spread, contrast)
%     sigmaG = (0.15/27)*spread + (0.25 - (0.15/27)*30 ) + ...
%         ( - (0.2/0.04)*contrast + (0.3 + (0.2/0.04)*0.01) + 0.4);
    sigmaG = (0.15/27)*spread + (0.25 - (0.15/27)*30 ) + ...
        ( - (0.2/0.04)*contrast + (0.3 + (0.2/0.04)*0.01) + 0.4);
    varGain = sigmaG.^2;
end

function [ll,xFit,pFit]  = computeLikelihood(outputLabels, predictors)
    b_coeff = glmfit(predictors(:), outputLabels(:), 'binomial'); % fit model
    probs   = glmval(b_coeff, predictors(:), 'logit'); % get predicted probabilities P(confReports=1∣Vc​=x)
    
    ll = sum(outputLabels(:).*log(probs(:)) + (1-outputLabels(:)).*log(1-probs(:)));

    xFit = linspace(min(predictors), max(predictors), 200);
    pFit = glmval(b_coeff, xFit, 'logit');

end

function losses = crossValidate(outputLabels, predictors)
    k = 5;
    cv = cvpartition(predictors(:), 'KFold', k);
    
    lossA = crossval(@(Xtrain,Ytrain,Xtest,Ytest) ...
        logLoss(glmfit(Xtrain, Ytrain, 'binomial'), Xtest, Ytest), ...
        predictors(:), outputLabels(:), 'partition', cv);
    
    losses = -lossA;
end

function L = logLoss(b, Xtest, Ytest)
    p = glmval(b, Xtest, 'logit');
    %eps = 1e-10;
    L = -sum(Ytest.*log(p) + (1-Ytest).*log(1-p)); % +eps
end