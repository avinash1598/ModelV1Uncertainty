clear all
close all

% 1. Estimate gain variance based on the defined formoula
% 2. Gain for each neuron comes from a distribution with mean gain variance
% defined by step 1 and variance of gain variance dependent on gain
% variance (linear relationship)
% 3. Higher peak firing rate, lower gain variability
% 4. Statistics of gain variance is preserved for a single neuron across
% different stimulus condition (i.e. correlated gain variances across 
% neurons for two stimulus condition)

% Maybe do bivariate without correlationt to firing rate

% 1. For each level obtain mean and variance of lognormal distribution
% 2. Convert to lognormal space
% 3. Build covariance matric for gains


nNeurons                = 200;  
tuningFnContrasts       = linspace(1e-4, 0.15, 49);
tuningFnSpreads         = linspace(1, 90, 50);
gains                   = zeros(numel(tuningFnContrasts), numel(tuningFnSpreads));
nrnGainMatrix           = zeros(nNeurons, numel(tuningFnContrasts), numel(tuningFnSpreads));

for cIdx=1:numel(tuningFnContrasts)
    for sIdx=1:numel(tuningFnSpreads)
        spread = tuningFnSpreads(sIdx);
        contrast = tuningFnContrasts(cIdx);
        gains(cIdx, sIdx) = getVarGain(spread, contrast);
        
        nrnGainMatrix(:,cIdx, sIdx) = getVarGainForAllNeuron(nNeurons, ...
            spread, contrast);
    end
end

figure
subplot(2, 2, 1)
histogram(nrnGainMatrix(:,1,1))

subplot(2, 2, 2)
histogram(nrnGainMatrix(:,49,50))

%%
function nrnsVarGain = getVarGainForAllNeuron(nNeurons, spread, contrast)
    varGain = getVarGain(spread, contrast);

    m = varGain;
    v = ( 0.4*m ).^2;
    sigma2 = log(1 + v./m.^2);
    sigma = sqrt(sigma2);
    mu = log(m) - 0.5*sigma2;
    
    nrnsVarGain = mu + sigma * randn(nNeurons,1);
    nrnsVarGain = exp(nrnsVarGain);
end

function varGain = getVarGain(spread, contrast)
    sigmaG = (0.15/27)*spread + (0.25 - (0.15/27)*30 ) + ...
        ( - (0.2/0.04)*contrast + (0.3 + (0.2/0.04)*0.01) + 0.4);
    varGain = sigmaG.^2;
end