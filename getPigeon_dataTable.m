function [dataTable_] = getPigeon_dataTable(options)
% function [dataTable_] = getPigeon_dataTable(options)
%
% Returns dataTable with the following columns:
%   1. subjectIndex
%   2. blockIndex
%   3. trialNumber
%   4. bound (pigeon position in OL)
%   5. boundLo (bottom of range)
%   6. boundHi (bottom of range)
%   7. RT (response time, including ndt, in number of steps)
%   8. DT (response time, not including ndt, in number of steps)
%   9. choice (0/1)
%   10. correct (correct=1,error=0)
%   11. coinCount
%   12. steps (cell)

arguments
    options.taskType = 'OL'    % 'OL', 'MX', or 'PD'
    options.ndt = 1            % non-decision time
    options.dataDir = ...
        fullfile('/Users', 'jigold', 'GoldWorks', 'Mirror_jigold', 'Manuscripts', '2023_Pigeon', 'Data')
    options.blocks = 'all';
end

% Set up the data table
dataTable_ = getPigeon_blankDataTable(0);

% Loop through the file list
files = dir(fullfile(options.dataDir, ['Pigeon_' options.taskType], 'prolificcsvs', '*.csv'));
for ff = 1:length(files)

    disp(ff)

    % Read data into Table
    tmpTable = readtable(fullfile(files(ff).folder, files(ff).name));

    % Use only good entries from selected blocks
    Lgood = isfinite(tmpTable.trial_number);
    if isnumeric(options.blocks) && sum(options.blocks) ~= sum(1:6)
        Lgood = Lgood & ismember(tmpTable.block_number, options.blocks);
    end

    % Get the columns we want
    if strcmp(options.taskType, 'OL') || strcmp(options.taskType, 'MX')
        tmpTable = tmpTable(Lgood, ...
            {'block_number', 'trial_number', 'pigeon_steps', 'choice', ...
            'correct' 'steps_count', 'step_std_val', 'coins_count'});
    else
        indices = find(Lgood);
        tmpTable = cat(2, ...
            tmpTable(indices, {'block_number', 'trial_number', 'predefined_bound', 'correct'}), ...
            tmpTable(indices+1, {'pigeon_steps', 'choice', 'coins_count'}));
    end

    % Add the data
    numTrials = size(tmpTable, 1);
    tableToAppend = getPigeon_blankDataTable(numTrials);
    tableToAppend.subjectIndex(1:end) = ff;
    tableToAppend.blockIndex(1:end) = tmpTable.block_number(1:end);
    tableToAppend.trialNumber(1:end) = tmpTable.trial_number(1:end);

    % Loop through the trials to fill in the blanks
    for tt = 1:numTrials

        % Get steps
        tableToAppend.steps{tt} = str2num(tmpTable.pigeon_steps{tt});

        % Get bounds, etc
        if strcmp(options.taskType, 'OL') || strcmp(options.taskType, 'MX')
            vals = getPigeon_bounds(tableToAppend.steps{tt}, options.ndt);
            tableToAppend.RT(tt) = length(tableToAppend.steps{tt});
            tableToAppend.bound(tt) = vals(1);
            tableToAppend.boundLo(tt) = vals(2);
            tableToAppend.boundHi(tt) = vals(3);
            tableToAppend.DT(tt) = vals(4);
        else
            tableToAppend.bounds(tt) = tmpTable.predefined_bound(tt);
            tableToAppend.rt(tt) = length(tableToAppend.steps(tt));
        end

        % Choice/outcome
        if ismember(tmpTable.choice{tt}, {'right', 'left'})
            tableToAppend.choice(tt) = double(strcmp(tmpTable.choice{tt}, 'right'));
            tableToAppend.correct(tt) = double(strcmp(tmpTable.choice{tt}, ...
                tmpTable.correct{tt}));
        end

        % Coin count
        tableToAppend.coinCount(tt) = tmpTable.coins_count(tt);
    end

    % add the data
    dataTable_ = cat(1, dataTable_, tableToAppend);
end