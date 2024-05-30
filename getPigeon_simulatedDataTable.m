function dataTable_ = getPigeon_simulatedDataTable(specs)
%
% Returns:
%   dataTable_ as in getPigeon_dataTable
%
% Arguments:
%   specs is either:
%       dataTable     ... table of human data to mimic
%       blockStructs  ... array of structs with (optional) substructs:
%           blockArgs
%           simArgs

% Defaults
blockDefaults = struct( ...
    'blockIndex',            nan, ... % if 1-6, use defaults
    'numSubjects',           60, ...
    'stepsPerBlock',         600, ...
    'P2Psteps',              1, ...
    'P2Pcoins',              0, ...
    'coinsGainedPerCorrect', 1, ...
    'coinsLostPerError',     0, ...
    'stepsLostPerError',     0);

simDefaults = struct( ...
    'generativeMean',        0.05, ...
    'generativeSTD',         0.15, ...
    'numTrials',             610, ...
    'maxStepsPerTrial',      50, ...
    'boundMean',             0, ...
    'boundSTD',              0, ...
    'boundSlope',            0, ...
    'boundMax',              0.75, ...
    'NDTMin',                1, ...
    'NDTmax',                3, ...
    'lapseRate',             0.05);

% Parse human data
if istable(specs)
    subjects = nonanunique(specs.subjectIndex);
    numSubjects = length(subjects);
    blocks = nonanunique(specs.blockIndex);
    numBlocks = length(blocks);
    Lgood = specs.trialNumber & specs.RT>=0; % non-zero RT
    bounds = nan(numSubjects,2);
    blockStructs = struct;
    for bb = 1:numBlocks
        for ss = 1:numSubjects
            Lsb = Lgood & specs.blockIndex==blocks(bb) & specs.subjectIndex==subjects(ss);
            bounds(ss,:) = [median(abs(specs.bound(Lsb)),'omitnan') ...
                std(abs(specs.bound(Lsb)),'omitnan')];
        end
        blockStructs(bb).blockArgs.blockIndex = bb;
        blockStructs(bb).simArgs.boundMean = bounds(:,1);
        blockStructs(bb).simArgs.boundSTD = bounds(:,2);
    end
else
    blockStructs = specs;
end

% Parse the arguments so we can estimate the table size
blockArgs = repmat(blockDefaults,1,length(blockStructs));
simArgs = repmat(simDefaults,1,length(blockStructs));
numTrials = 0;
for bb = 1:length(blockStructs)

    % fprintf('Simulation block %d\n', bb)

    % Get block arguments
    if isfield(blockStructs(bb), 'blockArgs')
        for ff = fieldnames(blockStructs(bb).blockArgs)'
            blockArgs(bb).(ff{:}) = blockStructs(bb).blockArgs.(ff{:});
        end
    end

    % Get simulation arguments
    if isfield(blockStructs(bb), 'simArgs')
        for ff = fieldnames(blockStructs(bb).simArgs)'
            simArgs(bb).(ff{:}) = blockStructs(bb).simArgs.(ff{:});
        end
    end

    % Check for blockIndex
    if blockArgs(bb).blockIndex>=1 && blockArgs(bb).blockIndex<=6
        switch blockArgs(bb).blockIndex
            case 1
                simArgs(bb).generativeMean = 0.05;
                simArgs(bb).generativeSTD = 0.15;
                blockArgs(bb).stepsPerBlock = 600;
                blockArgs(bb).P2Psteps = 1;
                blockArgs(bb).P2Pcoins = 0;
                blockArgs(bb).coinsGainedPerCorrect = 1;
                blockArgs(bb).coinsLostPerError = 0;
                blockArgs(bb).stepsLostPerError = 0;
            case 2
                simArgs(bb).generativeMean = 0.05;
                simArgs(bb).generativeSTD = 0.15;
                blockArgs(bb).stepsPerBlock = 600;
                blockArgs(bb).P2Psteps = 1;
                blockArgs(bb).P2Pcoins = 0;
                blockArgs(bb).coinsGainedPerCorrect = 1;
                blockArgs(bb).coinsLostPerError = 4;
                blockArgs(bb).stepsLostPerError = 0;
            case 3
                simArgs(bb).generativeMean = 0.05;
                simArgs(bb).generativeSTD = 0.15;
                blockArgs(bb).stepsPerBlock = 600;
                blockArgs(bb).P2Psteps = 1;
                blockArgs(bb).P2Pcoins = 0;
                blockArgs(bb).coinsGainedPerCorrect = 1;
                blockArgs(bb).coinsLostPerError = 0;
                blockArgs(bb).stepsLostPerError = 30;
            case 4
                simArgs(bb).generativeMean = 0.15;
                simArgs(bb).generativeSTD = 0.15;
                blockArgs(bb).stepsPerBlock = 600;
                blockArgs(bb).P2Psteps = 1;
                blockArgs(bb).P2Pcoins = 0;
                blockArgs(bb).coinsGainedPerCorrect = 1;
                blockArgs(bb).coinsLostPerError = 0;
                blockArgs(bb).stepsLostPerError = 0;
            case 5
                simArgs(bb).generativeMean = 0.15;
                simArgs(bb).generativeSTD = 0.15;
                blockArgs(bb).stepsPerBlock = 600;
                blockArgs(bb).P2Psteps = 1;
                blockArgs(bb).P2Pcoins = 0;
                blockArgs(bb).coinsGainedPerCorrect = 1;
                blockArgs(bb).coinsLostPerError = 4;
                blockArgs(bb).stepsLostPerError = 0;
            case 6
                simArgs(bb).generativeMean = 0.15;
                simArgs(bb).generativeSTD = 0.15;
                blockArgs(bb).stepsPerBlock = 600;
                blockArgs(bb).P2Psteps = 1;
                blockArgs(bb).P2Pcoins = 0;
                blockArgs(bb).coinsGainedPerCorrect = 1;
                blockArgs(bb).coinsLostPerError = 0;
                blockArgs(bb).stepsLostPerError = 30;
        end
    end


    % Update estimate of numTrials
    simArgs(bb).numTrials = blockArgs(bb).stepsPerBlock + 2;
    numTrials = numTrials + blockArgs(bb).numSubjects .* simArgs(bb).numTrials;
end

% Set up the data table, steps celery
%   Note that we're using the same fields as from getDataTable, with the
%   addition of "trueBound" -- bound is computed empirically, trueBound is,
%   duh, the real one
% Set up the data table -- preallocate for speed
dataTable_ = getPigeon_blankDataTable(numTrials, 'trueBound', "double");
tableIndex = 1;

% Now loop through the blocks and do the simulations
for bb = 1:length(blockStructs)

    % To convert struct to celery for args
    argNames = fieldnames(simArgs(bb));
    argVals = struct2cell(simArgs(bb));
    args = cell(1,length(argNames)*2);
    args(1:2:end) = argNames;

    stepCounts = nan(simArgs(bb).numTrials,1);
    coinCounts = nan(simArgs(bb).numTrials,1);
    NDT = min(2, simArgs(bb).NDTmax);

    % Do the simulation, per subject
    for ss = 1:blockArgs(bb).numSubjects

        % Make the arg list for the simulations
        for vv = 1:length(argVals)
            if length(argVals{vv}) == blockArgs(bb).numSubjects
                args{vv*2} = argVals{vv}(ss);
            else
                args{vv*2} = argVals{vv}(1);
            end
        end

        % Do the sim
        [choices, rts, DV, bounds, ~] = getPigeon_simulatedData(args{:});

        % Update the data Table
        % jig changed to match real data, where steps always include zero
        % start
        %         stepCounts(:) = cumsum(repmat(blockArgs(bb).P2Psteps,simArgs(bb).numTrials,1) + ...
        %             steps + blockArgs(bb).stepsLostPerError .* (choices==0));
        stepCounts(:) = cumsum(rts + blockArgs(bb).stepsLostPerError .* (choices==0));
        coinCounts(:) = cumsum(repmat(-blockArgs(bb).P2Pcoins,simArgs(bb).numTrials,1) - ...
            blockArgs(bb).coinsLostPerError .* (choices==0) + ...
            blockArgs(bb).coinsGainedPerCorrect .* (choices==1));
        trialCount = find(stepCounts>=blockArgs(bb).stepsPerBlock,1);

        % Add for final aborted trial
        if stepCounts > blockArgs(bb).stepsPerBlock
            choices(trialCount) = nan;
            coinCounts(trialCount) = coinCounts(trialCount-1);
        end

        % Make the table
        inds = tableIndex:tableIndex+trialCount-1;
        if inds(end) > height(dataTable_)
            dataTable_ = cat(1, dataTable_, ...
                getPigeon_blankDataTable(10000, 'trueBound', "double"));
        end

        dataTable_.subjectIndex(inds) = ss;
        dataTable_.blockIndex(inds) = bb;
        dataTable_.trialNumber(inds) = 1:trialCount;
        dataTable_.RT(inds) = rts(1:trialCount);
        dataTable_.trueBound(inds) = bounds(1:trialCount);
        dataTable_.choice(inds) = choices(1:trialCount);
        dataTable_.correct(inds) = choices(1:trialCount);
        dataTable_.coinCount(inds) = coinCounts(1:trialCount);

        % Compute lower/upper/bound + ndt
        rts = rts(1:trialCount);
        DV  = DV(1:trialCount,:);

        % Special case of one step
        Lrt = rts == 1;
        if any(Lrt)
            dataTable_.bound(inds(Lrt)) = DV(Lrt,1)./2;
            dataTable_.boundLo(inds(Lrt)) = 0;
            dataTable_.boundHi(inds(Lrt)) = DV(Lrt,1);
            dataTable_.DT(inds(Lrt)) = 1;
            dataTable_.steps(inds(Lrt)) = num2cell(DV(Lrt,1));
        end

        % General case
        for rr = 2:max(rts)
            Lrt = rts == rr;
            if any(Lrt)
                r_ndt = min(NDT, rr-2);
                dataTable_.bound(inds(Lrt)) = mean(DV(Lrt,rr-r_ndt-1+[0 1]),2);
                dataTable_.boundLo(inds(Lrt)) = DV(Lrt,rr-r_ndt-1);
                dataTable_.boundHi(inds(Lrt)) = DV(Lrt,rr-r_ndt);
                dataTable_.DT(inds(Lrt)) = rr-r_ndt;
                dataTable_.steps(inds(Lrt)) = mat2cell(DV(Lrt,1:rr), ones(sum(Lrt),1), rr);
            end
        end

        % update index
        tableIndex = inds(end)+1;
    end
end

% Remove fluff
dataTable_ = dataTable_(1:tableIndex-1,:);
