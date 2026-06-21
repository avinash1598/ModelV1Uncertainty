% ---------------------------------------------------------------------
% Function to decode stimulus orientation from spike data using a Bayesian approach.
% This function performs decoding based on the neuron's tuning curves and 
% the observed spike counts, estimating the most likely stimulus orientation.

% INPUTS:
% - spikes:            [nNeurons x nTimeBins] Binary matrix containing the spike responses for all neurons over time.
% - neuronsPrefOrientation: [1 x nNeurons] Preferred orientations of the neurons (in radians), determining their tuning.
% - params:            Structure containing several parameters related to the time bins, number of neurons, and stimulus duration.
%   - params.timeBins: [1 x nTimeBins] Time bin vector representing the time points for spike recording.
%   - params.nNeurons: Scalar representing the total number of neurons.
%   - params.stimDuration: Scalar representing the total duration of stimulus presentation (in seconds).
% - tuningParams:      Structure containing tuning curve parameters for the population of neurons (used for computing firing rates based on orientation).
%   - tuningParams: Includes various fields such as alpha, q, beta, etc., which control the neuron tuning properties.

% OUTPUTS:
% - thetaMLE: [1 x numIntervals] The maximum likelihood estimate (MLE) of the stimulus orientation for each time window.
% ---------------------------------------------------------------------
function [contrastMLE, spreadMLE, thetaMLE, pdfData, MLEs, metrics] = decodeModulatedPoissonSpikes(spikes, tuningFns, params, sigma_g2_all)
    % Decode assuming poisson spiking
    
    contrasts    = params.contrasts;
    spreads      = params.spreads;
    uniqStimOris = params.uniqStimOris;
    stimDuration = params.stimDuration;     % Stimulus duration (total time)
    
    %dOri = rad2deg( uniqStimOris(2) ) - rad2deg( uniqStimOris(1) );
    dOri = uniqStimOris(2) - uniqStimOris(1);
    dC   = contrasts(2) - contrasts(1);
    dS   = spreads(2) - spreads(1);

    % Decode orientation from aggrgate spike count (no need to consider small individual time windows)
    aggregateSpikeCounts = sum(spikes, 2);
    
    logPdf3D = zeros( numel(contrasts), numel(spreads), numel(uniqStimOris) );

    for cIdx = 1:numel(contrasts)
        for sIdx = 1:numel(spreads)

            key1 = sprintf('c_%g', contrasts(cIdx));
            key2 = sprintf('s_%g', spreads(sIdx));
            key1 = matlab.lang.makeValidName(key1);
            key2 = matlab.lang.makeValidName(key2);
            
            nrnsOriTuningFn = tuningFns.(key1).(key2); % For all neurons
            sigma_g2         = sigma_g2_all.(key1).(key2);
            sigma_g2         = sigma_g2'; % nNeurons x 1
            
            % compute PDF
            n = aggregateSpikeCounts;     % Spike count vector for all neurons in the current window
            f_theta = nrnsOriTuningFn;        % Firing rates for all neurons for each candidate orientation
            
            % Add a small epsilon to the firing rates to avoid numerical issues (log of 0)
            % TODO: make sure adding epsilon is okay
%             temp = sigma_g.^2 .* stimDuration .* f_theta + eps; 
%             log_pdf_p1 = log( gamma(n + 1./sigma_g.^2) ) - log(gamma(n + 1)) - log(gamma(1./sigma_g.^2));
%             log_pdf_p2 = n.*log(temp) - (n + 1./sigma_g.^2).*log(1 + temp);
%             log_pdf = sum(log_pdf_p1 + log_pdf_p2, 1, 'omitnan'); 
            
            temp = sigma_g2 .* stimDuration .* f_theta; % + eps (eps is problamatic)

            log_pdf_p1 = gammaln(n + 1./sigma_g2) ...
                       - gammaln(n + 1) ...
                       - gammaln(1./sigma_g2);
            
            log_pdf_p2 = n .* log(temp) ...
                       - (n + 1./sigma_g2) .* log1p(temp);  % more stable than log(1+temp)
            
            log_pdf = sum(log_pdf_p1 + log_pdf_p2, 1); % 'omitnan'
            
            if any(~isfinite(log_pdf(:)))  
                % true if NaN or Inf
                keyboard
            end
            
            logPdf3D(cIdx, sIdx, :) = log_pdf;

            %keyboard
        end
    end
    
    [~, idx] = max(logPdf3D(:));   % max value + linear index
    [x, y, z] = ind2sub(size(logPdf3D), idx);
    
    contrastMLE_v1 = contrasts(x);
    spreadMLE_v1   = spreads(y);
    thetaMLE_v1    = uniqStimOris(z);
    
    % pdfData = logPdf3D; % nContrast x nSpreads x nOris
    
    % convert log PDF to normalized PDF (use rad here)
    %pdf = exp(logPdf3D);
    %pdf = pdf./sum(pdf*dOri*dC*dS, 'all');
    
    %pdfOri      = squeeze(sum(sum(pdf,2),1)) * dC * dS;
    %pdfContrast = squeeze(sum(sum(pdf,3),2)) * dOri * dS;
    %pdfSpread   = squeeze(sum(sum(pdf,3),1)) * dOri * dC;
    
    %pdf = exp(logPdf3D);
    pdf = exp(logPdf3D - max(logPdf3D(:)));
    normFactor = trapz(contrasts, ...
                    trapz(spreads, ...
                    trapz(uniqStimOris, pdf, 3), 2), 1);
    pdf = pdf / normFactor;
    
    pdfOri      = squeeze(trapz(contrasts, trapz(spreads, pdf, 2), 1));
    pdfContrast = squeeze(trapz(spreads, trapz(uniqStimOris, pdf, 3), 2));
    pdfSpread   = squeeze(trapz(contrasts, trapz(uniqStimOris, pdf, 3), 1));
    
    % Keep PDFs and use different metrics to calculate std and mean
    %pdfData.logPdf3D    = logPdf3D;
    %pdfData.pdf3D       = pdf;
    pdfData.pdfOri      = pdfOri(:);
    pdfData.pdfContrast = pdfContrast(:);
    pdfData.pdfSpread   = pdfSpread(:);
    
    [~, idx] = max(pdfOri(:));
    thetaMLE_v2 = uniqStimOris(idx);
    
    [~, idx] = max(pdfContrast(:));
    contrastMLE_v2 = contrasts(idx);
    
    [~, idx] = max(pdfSpread(:));
    spreadMLE_v2 = spreads(idx);
    
    MLEs.thetaMLE_v1    = thetaMLE_v1; % joint pdf
    MLEs.contrastMLE_v1 = contrastMLE_v1; % joint pdf
    MLEs.spreadMLE_v1   = spreadMLE_v1; % joint pdf
    % Use marginally decoded values instead
    MLEs.thetaMLE_v2    = thetaMLE_v2; % marginal pdf
    MLEs.contrastMLE_v2 = contrastMLE_v2; % marginal pdf
    MLEs.spreadMLE_v2   = spreadMLE_v2; % marginal pdf
    
    % Uncertainty computation
    theta = uniqStimOris(:);   % in radians
    w = pdfOri(:) ; %* dOri;      % weights
    C = trapz(theta, w .* cos(2*theta));
    S = trapz(theta, w .* sin(2*theta));
    %C = sum(w .* cos(2*theta));
    %S = sum(w .* sin(2*theta));

    mu    = 0.5 * atan2(S, C); mu = mod(mu, pi);
    R     = sqrt(C^2 + S^2); % V = 1 - R; % circular variance
    sigma = sqrt(-2 * log(R)) / 2; % circular standard deviation
    
    mu2    = sum( uniqStimOris(:).*pdfOri(:)*dOri );
    sigma2 = sqrt( sum( (uniqStimOris(:) - mu2).^2.*pdfOri(:)*dOri ) );
    
    metrics.mu2 = mu2;
    metrics.sigma2 = sigma2;
    metrics.mu = mu;
    metrics.sigma = sigma;
    
    % Return the marginally decoded value
    contrastMLE = contrastMLE_v2;
    spreadMLE   = spreadMLE_v2;
    thetaMLE    = thetaMLE_v2;
    
    %keyboard
end