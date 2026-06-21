clear all
% close all

% Stim profile
% No direction
% srtimulus energy is only from 0 to 180 degrees

contrasts = [0.01 1];
dispersions = [3 30];
oriSpace = linspace(-pi, pi, 501); % linspace(0, pi, 501)

figure

for c = 1:numel(contrasts)
    for d = 1:numel(dispersions)
        stimParams.contrastLevel = contrasts(c);
        stimParams.spreadLevel   = deg2rad(dispersions(d));
        stimParams.stimOri       = deg2rad(70);

        sp = generateStimProfile(oriSpace, stimParams);
        dtheta = oriSpace(2) - oriSpace(1);
        energy = sum( (sp.^2) .* dtheta );
        
        subplot(2, 2, c)
        hold on
        plot(rad2deg(oriSpace), sp, DisplayName=sprintf("D=%d", dispersions(d)))
        title( sprintf("Total contrast energy =%.2f", contrasts(c)) )
        xlabel("orientation")
        ylabel("Contrast")
        hold off
        legend

    end
end

