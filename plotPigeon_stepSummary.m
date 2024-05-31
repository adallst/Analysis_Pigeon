function stepData = plotPigeon_stepSummary(positions, options)

arguments
    positions               cell
    options.choices         double = []
    options.LgoodTrials     logical = []
    options.groups          double = []      
    options.maxSteps        double = 50   % Maximum steps per trial
    options.numSteps        double = 10   % Number steps to show
    options.axs             double = gca
    options.generativeMean  double = 0.05
    options.titletext           string = "Simulation"
end

% Make steps matrix
numTrials = size(positions, 1);
stepMatrix = nan(numTrials, options.maxSteps);
for xx = 1:numTrials
    stepMatrix(xx,1:length(positions{xx})-1) = flip(diff(positions{xx}));
end

% Possibly sign by choice
if ~isempty(options.choices)
    stepMatrix = stepMatrix.*repmat(options.choices,1,options.maxSteps);
end

% Update selection array
if isempty(options.LgoodTrials)
    options.LgoodTrials = true(numTrials, 1);
end

% Compute mean step size per step re: end
if isempty(options.groups)
    options.groups = ones(numTrials, 1);
end
uniqueGroups = nonanunique(options.groups);
numGroups = length(uniqueGroups);
stepData = nan(numGroups, options.numSteps);
for gg = 1:numGroups
    Lg = options.LgoodTrials & options.groups == uniqueGroups(gg);
    if any(Lg)
        stepData(gg,:) = mean(stepMatrix(Lg,1:options.numSteps),'omitnan');
    end
end

% Plot it
xax = (-1:-1:-options.numSteps)';
axes(options.axs); cla reset; hold on;
plot([-options.numSteps -1], [0 0], 'k:')
plot([-options.numSteps -1], [options.generativeMean options.generativeMean], 'b--')
if numGroups > 2
    plot(xax, stepData', '-', 'Color', 0.8.*ones(3,1))
    plot(xax, mean(stepData,'omitnan'), 'r-', 'LineWidth', 3)
else
    plot(xax, stepData, 'r-', 'LineWidth', 3)
end
axis([-options.numSteps -1 -0.05 0.16])
xlabel('Past steps re: choice')
ylabel('Step size')
title(options.titletext)
