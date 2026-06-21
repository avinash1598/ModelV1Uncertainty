% Keep stimulus energy constant for one contrast.

function Sp = generateStimProfile(oriSpace, stimParams)
    % tuningFnOrientations : in radians
    % A direction parameter maybe???

    contrastLevel = stimParams.contrastLevel';
    spreadLevel   = stimParams.spreadLevel';
    stimOri       = stimParams.stimOri';
    
    % Keep the contrast energy constant
    netContrastEnergy = contrastLevel.^2;
    dtheta = oriSpace(2) - oriSpace(1);

    % scale     = contrastLevel;
    angleDiff = getAcuteAngleDiff(oriSpace, stimOri);
    sigma     = spreadLevel;
    
    % Sp = scale.*exp( -angleDiff.^2 ./ (2*sigma.^2) ); % ntrials x nOris
    Sp = exp( -angleDiff.^2 ./ (2*sigma.^2) ); % ntrials x nOris

    % Compute contrast energy 
    energy = sum( (Sp.^2) .* dtheta );

    % Normalized stimulus profile
    Sp = Sp .* sqrt( ( netContrastEnergy / energy ) ); % This does not change the dispersion

    % keyboard
end

function d = getAcuteAngleDiff(ori1, ori2)
    d = angle(exp(1i * 2*(ori1 - ori2))) / 2;
    % d = rad2deg(d);
end
