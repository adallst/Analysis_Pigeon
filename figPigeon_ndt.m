function figPigeon_ndt(dataTable, simDataTable, generativeMean, num)
% function figPigeon_ndt(dataTable, simDataTable, generativeMean, num)
%
% Figure: non-decision times
%
arguments
    dataTable    
    simDataTable
    generativeMean = 0.5
    num = 1.1
end

%% Set up figure
%
wid     = 8.5; % total width
cols    = {1,1,1};
hts     = [3.5 3.5 3.5];
[axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
set(axs,'Units','normalized');

%% 1. Simulate basic DDM, no ndt
%
blockStruct.blockArgs.stepsPerBlock = 10000;
blockStruct.simArgs.boundMean = 0.4;
blockStruct.simArgs.NDTMin = 0;
blockStruct.simArgs.NDTMax = 0;
sdt = getPigeon_simulatedDataTable(blockStruct);
plotPigeon_stepSummary(sdt.steps, ...
    'axs',              axs(1), ...
    'LgoodTrials',      sdt.choice==1, ...
    'generativeMean',   generativeMean);

%% 2. Use data-matched simulations
%
Lgood = simDataTable.trialNumber>5 & simDataTable.correct==1; % Correct only!
stepData = plotPigeon_stepSummary(simDataTable.steps, ...
    'axs',              axs(2), ...
    'groups',           simDataTable.subjectIndex, ...
    'choices',          simDataTable.choice.*2-1, ...
    'LgoodTrials',      Lgood, ... 
    'generativeMean',   generativeMean, ...
    'titletext',        "Data-matched simulations");

%% Data
%

% Get the standard data matrix
Lgood = dataTable.trialNumber>5 & dataTable.correct==1; % Correct only!
stepData = plotPigeon_stepSummary(dataTable.steps, ...
    'axs',              axs(3), ...
    'groups',           dataTable.subjectIndex, ...
    'choices',          dataTable.choice.*2-1, ...
    'LgoodTrials',      Lgood, ... 
    'generativeMean',   generativeMean, ...
    'titletext',        "Data");

% % % plot step size 1&2 back vs 0 back
% % wh = 0.99.*ones(3,1);
% % subplot(4,1,4); cla reset; hold on;
% % plot([-0.02 0.2], [-0.02 0.2], 'k:')
% % plot(stepData(:,2)-stepData(:,1), stepData(:,3)-stepData(:,1), 'ko', ...
% %     'MarkerFaceColor', wh)
% plot(stepData(:,[1 1])', stepData(:,[2 3])', 'k-');
% plot(stepData(:,1), stepData(:,2), 'ro', 'MarkerFaceColor', 'r');
% plot(stepData(:,1), stepData(:,3), 'bo', 'MarkerFaceColor', 'b');
