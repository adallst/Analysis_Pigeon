%% Collect data Tables

% Real data
dataTable = getPigeon_dataTable( ...
    'taskType',     'OL', ... % 'OL', 'MX', or 'PD' (not tested with PD)
    'blocks',       1:3);

% Simulated data
% 
% Includes random non-decision time of 1-3 steps
% Uses median/STD bounds from each subject to do simulated subjects
simDataTable = getPigeon_simulatedDataTable(dataTable);

%% Figure 2: Performance Summary
%
% Simulations are a little slow, change arg to 10 for fast rendering
figPigeon_performanceSummary(dataTable, 2.0, 100)
figPigeon_performanceSummary(simDataTable, 2.1, 100)

%% Figure 3: Accuracy vs RT
%
figPigeon_accuracyVsRT(dataTable, 3.0)
figPigeon_accuracyVsRT(simDataTable, 3.1)

%% Figure 4: Congruence betweeen choice and pigeon position
% 
% shown for different non-decision times (rows)
figPigeon_choiceCongruence(dataTable, 4.0)
figPigeon_choiceCongruence(simDataTable, 4.1)

%% Figure 5: Normalized bound vs RT
% 
% Normalized to mean bound per subject
% Shown for different non-decision times (rows)
figPigeon_boundByTime(dataTable, 5.0)
figPigeon_boundByTime(simDataTable, 5.1)

%% Figure 6: Compare bounds early/mid/late
% 
figPigeon_boundTertiles(dataTable, 6.0)
figPigeon_boundTertiles(simDataTable, 6.1)

%% Figure 7: Bound fits
%
figPigeon_boundFits(dataTable, 7.0)
figPigeon_boundFits(simDataTable, 7.1)

%% Figure 8: Learning
%
figPigeon_blockLearning(dataTable, 8.0)
figPigeon_blockLearning(simDataTable, 8.1)

%% Figure 9: Regression to mean
%
figPigeon_regressionToMean(dataTable, 9.0)
figPigeon_regressionToMean(simDataTable, 9.1)

%% Figure 1.S1: Computing NDT
%
figPigeon_ndt(dataTable, simDataTable, 1.1)

