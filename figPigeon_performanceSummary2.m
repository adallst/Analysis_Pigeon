function figPigeon_performanceSummary2(dataTable, block_names_publish, num, numSimSubjects)
% function figPigeon_performanceSummary(dataTable, num)
%
% Figure: basic performance
%
% top row: RT vs accuracy (% correct), per subject
% bottom row: coins/step vs median bound, per subject

arguments
    dataTable 
    block_names_publish
    num = 2
    numSimSubjects = 1000
    
end

%% Set up figure
wid     = 17.6./3; % total width
cols    = {1,1};
hts     = [3.5 3.5];
[axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);
pData = nan(numSubjects,6,numBlocks); % 2: pct correct, 25/50/75% RT, median bound, total reward

% Loop through each subject
fprintf('figPigeon_performanceSummary: Collecting real data\n')
Lcor = dataTable.correct==1;
Lg   = dataTable.RT>=0;
for ss = 1:numSubjects
    for bb = 1:numBlocks
        % Get subject- and block-specific data
        Lsb = Lg & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);

        % Pct correct
        pData(ss,1,bb) = sum(Lsb&Lcor)/sum(Lsb).*100;

        % RT median, IQR
        pData(ss,2:4,bb) = prctile(dataTable.RT(Lsb),[25 50 75]);

        % Median bound
        pData(ss,5,bb) = median(abs(dataTable.bound(Lsb)),'omitnan');

        % Total reward
        pData(ss,6,bb) = dataTable.coinCount(find(Lsb,1,'last'));
    end
end

%% Simulations
bounds = 0:0.01:0.75;
numBounds = length(bounds);
blockStructs = struct();
blockStructs.blockArgs.numSubjects = numSimSubjects;
rrs = nan(numBounds,3,numBlocks);
for bb = 2%1:numBlocks
    
    fprintf('figPigeon_performanceSummary: Collecting simulated data, block %d\n', bb)

    % Simulation parameters, no NDT
    blockStructs.blockArgs.blockIndex = bb;
    blockStructs.simArgs.NDTMin = 0;
    blockStructs.simArgs.NDTmax = 0;
   
    for dd = 1:numBounds

        % Set the bound
        blockStructs.simArgs.boundMean = bounds(dd);

        % Run the simulation
        sdt = getPigeon_simulatedDataTable(blockStructs);

        %fprintf('Computing optimal RR per bound, block %d, bound=%.2f, mean rt=%.2f, pct error=%.2f\n', ...
        %    bb, bounds(dd), mean(sdt.RT,'omitnan'), sum(sdt.correct==0)./sum(isfinite(sdt.correct)).*100)

        % Compute RR percentiles
        totalCoinsPerSimSubject = nan(numSimSubjects,1);
        for ss = 1:numSimSubjects
            lastTrialIndex = find(sdt.subjectIndex==ss,1,'last');
            totalCoinsPerSimSubject(ss) = sdt.coinCount(lastTrialIndex);
        end
        rrs(dd,:,bb) = prctile(totalCoinsPerSimSubject./600,[2.5 50 97.5]);
    end
end

%% Plotz
% Per block
%STEPS_PER_BLOCK = 600;
gry = 0.9.*ones(1,3);
wht = 0.99.*ones(1,3);
for bb = 2%1:numBlocks

    % RT vs pct correct
    axes(axs(1)); cla reset; hold on;
    plot(pData(:,[1 1],bb)', pData(:,[2 4],bb)', 'k-');
    plot(pData(:,1,bb), pData(:,3,bb), 'ko', 'MarkerFaceColor', wht);
    plot([40 100], median(pData(:,3,bb),'omitnan').*[1 1], 'b:')
    axis([40 100 0 18])
%     if bb == 1
        xlabel('Pct correct')
        ylabel('RT (steps)')
%     end
    title(block_names_publish(bb)) 

    % coins/step
    axes(axs(2)); cla reset; hold on;
    h=patch([bounds flip(bounds)], [rrs(:,1,bb)' flip(rrs(:,3,bb)')], gry);
    % sd = std(rrs(:,:,bb),'omitnan');
    % h=patch([bounds flip(bounds)], [mn-sd flip(mn+sd)], gr);
    set(h, 'LineStyle', 'none')
    plot(bounds, rrs(:,2,bb), 'k-', 'LineWidth', 2)

    plot(pData(:,5,bb), pData(:,6,bb)./600, 'ko', 'MarkerFaceColor', wht);
    plot(median(pData(:,5,bb)).*[1 1], [-1 1], 'b:')
    axis([0 0.75 -0.4 0.4])
%     if bb == 2
        xlabel('Median bound')
        ylabel('Coins/step')
%     end
end
