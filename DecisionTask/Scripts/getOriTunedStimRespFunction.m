function [tuningFnData] = getOriTunedStimRespFunction(neuronsPrefOrientations, params, stimParams)
% Note: there is no source of randomness in this function. It's compltely
% deterministic
    % ORIENTATIONTUNEDFIRINGRATE Computes the firing rates of neurons based on their 
    % preferred orientation using a tuning function model.
    %
    % INPUTS:
    % orientations - (nTrials x 1 vector) Array of stimulus orientations (in radians).
    %                These are the orientations that the subject is exposed to.
    %
    % neuronsPrefOrientations - (1 x nNeurons vector) Array of preferred orientations 
    %                           for each neuron (in radians). Each neuron is tuned to 
    %                           a specific orientation. 
    %
    % params - (struct) Structure containing the tuning function parameters:
    %   q    : Non-linearity exponent applied to the response of the neuron.
    %   eps1 : Baseline (offset) firing rate for each neuron.
    %   beta : Controls the range and scaling of the firing rates.
    %
    % stimParams - 
    %   contains trial level contrast, spread, and orientation information
    %   contrast
    %   dispersion
    % 
    % OUTPUT:
    % FR - (nTrials x nNeurons matrix) Firing rates of neurons in response to the given 
    %      stimulus orientations. Each row represents a trial (stimulus orientation), 
    %      and each column represents a neuron.
    %
    % DESCRIPTION:
    % The function computes the orientation-tuned firing rate of each neuron. This is
    % done by first calculating the unnormalized tuned response of each neuron to the 
    % given orientations, followed by normalization, and applying a power law 
    % non-linearity to adjust for firing rates.
    
    % All angles in radians

    % Orientation space over which tuning fn needs to be computed
    oriSpace = linspace(-pi, pi, 501);
    
    % Step 1: Get unnormalized tuning response of neurons based on their preferred orientation.
    unnormalizedResp = getUnormalizedTunedResponse( ...
        oriSpace, neuronsPrefOrientations, params); % nNeurons x nOris

    % Step2: Stim esponse f(S)
    stimProfile = generateStimProfile(oriSpace, stimParams); % nTrials x nOris
    stimfltResp = sum( unnormalizedResp.*stimProfile, 2 ); % Filter stimulus response 
    
    % Step 3: Compute the final firing rate using power law non-linearity.
    % This models how the response of the neuron changes based on non-linear scaling.
    q       = params.q;     % Transduction non-linearity exponent
    e1      = params.e1;    % Baseline firing rate (offset)
    e2      = params.e2;    % Baseline firing rate (stim dependent)
    gam     = params.gam;   % Response amplitude
    beta    = params.beta;  % Scaling factor
    % sigma_n = params.normalizationNoiseParam; % sigma_n for normalization noise
    
    % Is this the only randomness here?
    % normalizationNoise = 0 + sigma_n*randn(size(neuronsPrefOrientations));
    % assert( sum(normalizationNoise >= ( sum( stimfltResp ) ) ) == 0 );
    
    % Firing rates of this stimulus
    FR = e1 + gam.*( (e2 + stimfltResp)./( beta + sum( stimfltResp ) ) ).^q; % normalizationNoise
    % FR = e1 + gam.*( (e2 + stimfltResp)./(beta + mean( stimfltResp ) ) ).^q;
    
    tuningFnData.oriSpace             = oriSpace;
    tuningFnData.unnormalizedResp     = unnormalizedResp;
    tuningFnData.stimfltResp          = stimfltResp;
    tuningFnData.stimProfile          = stimProfile;
    tuningFnData.FR                   = FR;
    
end

function [FR] = getUnormalizedTunedResponse(orientations, neuronsPrefOrientations, params)
    % GETUNORMALIZEDTUNEDRESPONSE Computes the raw (unnormalized) tuning function response
    % of neurons based on the stimulus orientations and their preferred orientations.
    %
    % INPUTS:
    % orientations - (nTrials x 1 vector) Array of stimulus orientations (in radians).
    %
    % neuronsPrefOrientations - (1 x nNeurons vector) Preferred orientations for each neuron (in radians).
    %
    % params - (struct) Structure containing the tuning function parameters:
    %   d    : Direction selectivity of neurons. (Set to zero in this case, no direction selectivity.)
    %   alpha: Aspect ratio. Controls the sharpness of the tuning curve.
    %   b    : Sharpness of the tuning curve (power of exponent).
    %   q    : Non-linearity exponent applied after tuning curve.
    %   w    : Weight of untuned filter amplitude.
    %   UNTUNED_FILTER_AMPL : Amplitude of untuned filter response (set to zero in this case).
    %
    % OUTPUT:
    % FR - (nTrials x nNeurons matrix) Unnormalized firing rate responses for each trial (orientation)
    %      and each neuron (based on its preferred orientation).
    %
    % DESCRIPTION:
    % This function calculates the raw, unnormalized firing rate of neurons based on their 
    % orientation tuning curves. It incorporates factors like direction selectivity, sharpness 
    % (alpha and b), and untuned responses (if applicable). The result is an unrectified response,
    % which will later be normalized and modified by non-linearity.
    
    % Extract parameters from the input struct
    d             = params.d;                    % Direction selectivity
    alpha         = params.alpha;                % Aspect ratio, controls the sharpness of the curve.
    b             = params.b;                    % Controls the sharpness of the tuning curve.
    w             = params.w;                    % Weight of untuned filter
    untunedFltAmp = params.UNTUNED_FILTER_AMPL;  % Untuned filter amplitude.
    
    % Step 1: Compute the cosine of the angular difference between stimulus orientation and 
    % the neuron's preferred orientation.
    t1 = orientations - neuronsPrefOrientations;  % Angular difference (nTrials x nNeurons)

    % Step 2: Apply direction selectivity (if non-zero, modifies the response).
    t2 = 1 + 0.5 * (sign(cos(t1)) - 1) .* d;  % Direction-selective factor
    
    % Step 3: Compute the tuning curve response.
    % This uses a combination of cosine tuning and exponential factors to simulate neuron selectivity.
    t3 = ( cos(t1) .* exp( -0.5 * (cos(t1).^2) .* (1 - alpha.^2) ) ).^b;  % Tuning response
    
    % Step 4: Compute the tuned and untuned filter responses.
    rTuned = t2 .* t3;                % Tuned linear filter response
    rUntuned = w * untunedFltAmp;     % Untuned filter response (**zero** in this case)
    
    % Step 5: Combine the tuned and untuned responses.
    combinedChOp = rTuned - rUntuned;  % Overall response
    
    % Step 6: Apply rectifier (set negative values to zero).
    combinedChOp(combinedChOp < 0) = 0;
    
    % Step 7: Return the final unnormalized firing rate response.
    FR = combinedChOp; % nNeurons x nOris
    
    % Step 8: Normalize between 0 and 1 (not sure if this is fully right though) 
    FR = ( FR - min(FR, [], 2) ) ./ ( max(FR, [], 2) - min(FR, [], 2) );
    
end
