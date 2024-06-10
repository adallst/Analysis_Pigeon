function figPigeon_regressionToMean2(dataTable, block_names_publish, fdat, num)
% function figPigeon_regressionToMean(dataTable, num)
%
% Figure: regression to mean
%
arguments
    dataTable   
    block_names_publish
    fdat
    num = 2
end

%% Set up figure
wid     = 17.6; % total width
cols    = {3};
hts     = 3.0;
% [axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
% set(axs,'Units','normalized');
figure
tiledlayout(3,1)

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);
rData = nan(numSubjects,3,numBlocks,2); % per sub per iter: regression slope, intercept, p-value
rmhead = [0 20]; %removing the first so many trials comparison to all trials

% For shuffle text
numShuffles = 1000;
shuffleData = nan(numShuffles,2);

for ct =1:length(rmhead)
% Loop through each subject
Lg = dataTable.trialNumber > rmhead(ct) & dataTable.RT>=0;
for bb = 2%1:numBlocks
    fprintf('figPigeon_regressionToMean: Collecting shuffled data, block %d\n', bb)
    for ss = 1:numSubjects

        % Get subject- and block-specific bounds
        Lsb = Lg & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
        bounds = abs(dataTable.bound(Lsb));
%         bounds = bounds(round(length(bounds)/4):end);

        % Get the real linear fit
        numTrials = length(bounds);
        X = cat(2, ones(numTrials,1), bounds);
        Y = diff(bounds);
        rData(ss,1:2,bb,ct) = X(1:end-1,:)\Y;
       
        % Do the shuffles
        for ff = 1:numShuffles
            X(:,2) = bounds(randperm(numTrials));
            Y = diff(X(:,2));
            shuffleData(ff,:) = X(1:end-1,:)\Y;
        end

        rData(ss,3,bb,ct) = sum(shuffleData(:,2)<rData(ss,2,bb,ct))./numShuffles;
    end
end

%% Plotz
% Per block
wht = 0.99.*ones(1,3);
for bb = 2%1:numBlocks

    % zbound median, IQR
%     axes(axs(bb)); cla reset; hold on;
%     plot(rData(:,1,bb), rData(:,2,bb), 'ko', 'MarkerFaceColor', wht);
%     Lsig = rData(:,3,bb)<0.005 | rData(:,3,bb)>0.995;
%     plot(rData(Lsig,1,bb), rData(Lsig,2,bb), 'ko', 'MarkerFaceColor', 'k');
%     plot([0 0.75], [-1 -1], 'k:');
%     axis([0 0.75 -2 1])
%     if bb == 1
%         xlabel('Bound')
%         ylabel('Delta bound')
%     end
    
     nexttile; cla reset; hold on;
    plot(rData(:,1,bb,ct), rData(:,2,bb,ct), 'ko', 'MarkerFaceColor', wht);
    Lsig = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
    plot(rData(Lsig,1,bb,ct), rData(Lsig,2,bb,ct), 'ko', 'MarkerFaceColor', 'k');
    plot([0 0.75], [-1 -1], 'k:');
    axis([0 0.75 -2 1])
    if bb == 2
        xlabel('Intercept')
        ylabel({['Remove first ' num2str(rmhead(ct)) ' trials'],'Slope'})
    end
    if ct == 1
        title(block_names_publish(bb))
    end
end
end
%% 

% Lg = dataTable.trialNumber > rmhead(ct) & dataTable.RT>=0;
for bb = 2%1:numBlocks
    fprintf('figPigeon_regressionToMean: Collecting shuffled data, block %d\n', bb)
    for ss = 1:numSubjects

        % Get subject- and block-specific bounds
        Lsb = dataTable.trialNumber > fdat(ss,3,bb) & dataTable.RT>=0 & ...
            dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
        bounds = abs(dataTable.bound(Lsb));
%         bounds = bounds(round(length(bounds)/4):end);
        
        if ~isempty(bounds)&length(bounds)>1
        % Get the real linear fit
        numTrials = length(bounds);
        X = cat(2, ones(numTrials,1), bounds);
        Y = diff(bounds);
        rData(ss,1:2,bb,ct) = X(1:end-1,:)\Y;
       
        % Do the shuffles
        for ff = 1:numShuffles
            X(:,2) = bounds(randperm(numTrials));
            Y = diff(X(:,2));
            shuffleData(ff,:) = X(1:end-1,:)\Y;
        end

        rData(ss,3,bb,ct) = sum(shuffleData(:,2)<rData(ss,2,bb,ct))./numShuffles;
        end
    end
end

%% Plotz
% Per block
wht = 0.99.*ones(1,3);
for bb = 2%1:numBlocks

    % zbound median, IQR
%     axes(axs(bb)); cla reset; hold on;
%     plot(rData(:,1,bb), rData(:,2,bb), 'ko', 'MarkerFaceColor', wht);
%     Lsig = rData(:,3,bb)<0.005 | rData(:,3,bb)>0.995;
%     plot(rData(Lsig,1,bb), rData(Lsig,2,bb), 'ko', 'MarkerFaceColor', 'k');
%     plot([0 0.75], [-1 -1], 'k:');
%     axis([0 0.75 -2 1])
%     if bb == 1
%         xlabel('Bound')
%         ylabel('Delta bound')
%     end
    
     nexttile; cla reset; hold on;
    plot(rData(:,1,bb,ct), rData(:,2,bb,ct), 'ko', 'MarkerFaceColor', wht);
    Lsig = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
    plot(rData(Lsig,1,bb,ct), rData(Lsig,2,bb,ct), 'ko', 'MarkerFaceColor', 'k');
    plot([0 0.75], [-1 -1], 'k:');
    axis([0 0.75 -2 1])
    if bb == 2
        xlabel('Intercept')
        ylabel({['Remove first tau trials'],'Slope'})
    end
end

end
