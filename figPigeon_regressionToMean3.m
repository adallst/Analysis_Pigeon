function figPigeon_regressionToMean3(dataTable, block_names_publish, fdat, num)
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
tiledlayout(2,3)
EXAMPLE_SUBJECT=12; %find(sum(fitsdata2(:,3,:)<20,3)==3)
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
for bb = 1:numBlocks
    fprintf('figPigeon_regressionToMean: Collecting shuffled data, block %d\n', bb)
    for ss = 1:numSubjects

        % Get subject- and block-specific bounds
        if ct ==1 %from rmhead
            Lsb = Lg & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
            bounds = abs(dataTable.bound(Lsb));
        elseif ct ==2 %from tau fit
           % Get subject- and block-specific bounds
            Lsb = dataTable.trialNumber > fdat(ss,3,bb) & dataTable.RT>=0 & ...
                dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb); 
            bounds = abs(dataTable.bound(Lsb));
        end
%         bounds = bounds(round(length(bounds)/4):end);
        if sum(Lsb)>2
        % Get the real linear fit
        numTrials = length(bounds);
        X = cat(2, ones(numTrials,1), bounds);
        Y = diff(bounds);
        rData(ss,1:2,bb,ct) = X(1:end-1,:)\Y;
        
        if ss==EXAMPLE_SUBJECT&ct==2
            egX{bb}=X(1:end-1,2);
            egY{bb}=Y;
        end
       
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
end


%% Plotz
% Per block
gr = 0.9.*ones(1,3);
wt = 0.99.*ones(1,3);
for bb = 1:numBlocks

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
    
%      nexttile; cla reset; hold on;
%     plot(rData(:,1,bb,ct), rData(:,2,bb,ct), 'ko', 'MarkerFaceColor', wht);
%     Lsig = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
%     plot(rData(Lsig,1,bb,ct), rData(Lsig,2,bb,ct), 'ko', 'MarkerFaceColor', 'k');
%     plot([0 0.75], [-1 -1], 'k:');
%     axis([0 0.75 -2 1])
%     if bb == 2
%         xlabel('Intercept')
%         ylabel({['Remove first tau trials'],'Slope'})
%     end
    
%     nexttile; cla reset;
%     ct =1;
%     bins = -1.5:0.125:0.1;
%     Lsig1 = [];
%     Lsig1 = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
%     ct=2;
%     Lsig2 = [];
%     Lsig2 = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
%     scatterhistogram(rData(:,2,bb,1), rData(:,2,bb,2),"GroupData",Lsig1+(Lsig2.*2))
%     xlabel('Slope w/ all trials')
%     ylabel('Slope w/o tau trials')
%     title(block_names_publish(bb))
%     xlim([-2 1])
%     ylim([-2 1])
    
%     nexttile;hold on
%     ct =1;
%     bins = -1.5:0.125:0.1;
%     Lsig1 = [];
%     Lsig1 = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
%     ct=2;
%     Lsig2 = [];
%     Lsig2 = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
%     plot(rData(:,2,bb,1), rData(:,2,bb,2), 'ko', 'MarkerFaceColor', wt);
%     plot(rData(Lsig1,2,bb,1), rData(Lsig1,2,bb,2), 'ko', 'MarkerFaceColor', gr);
%     plot(rData(Lsig2,2,bb,1), rData(Lsig2,2,bb,2), 'ko', 'MarkerFaceColor', 'k');
%     xlabel('Slope w/ all trials')
%     ylabel('Slope w/o tau trials')
%     title(block_names_publish(bb))
%     xlim([-2 1])
%     ylim([-2 1])

    nexttile(bb);hold on
    yline(0,'k:')
    plot(0:0.025:0.75,0.75:-0.05:-0.75,'k:')
    plot(egX{bb},egY{bb},'k.')
    plot(0:0.025:0.75,[0:0.025:0.75].*rData(EXAMPLE_SUBJECT,2,bb,ct)+rData(EXAMPLE_SUBJECT,1,bb,ct),'k')
    xlabel('Previous bound')
    ylabel('Change in bound')
    title(block_names_publish(bb))

    bins = -1.5:0.125:0.1;
    nexttile(bb+numBlocks);hold on
    ct =2;
%     plot(rData(:,1,bb,ct), rData(:,2,bb,ct), 'ko', 'MarkerFaceColor',wt);
    Lsig = [];
    Lsig = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
%     plot(rData(Lsig,1,bb,ct), rData(Lsig,2,bb,ct), 'ko','MarkerFaceColor', 'k'); plot([0 0.75], [-1 -1], 'k:'); axis([0 0.75 -2 1])
    h1 = histcounts(rData(Lsig,2,bb,ct),'BinEdges',bins);
    h2 = histcounts(rData(~Lsig,2,bb,ct),'BinEdges',bins);
    bar(bins(2:end),[h1;h2]'./numSubjects,'stacked')
    newcolors = [0 0 0; 1 1 1];
    colororder(newcolors)
%     if bb == 3
%         xlabel('Intercept') %ylabel('Slope')
        xlabel('Slope')
%     end
    if bb ==1
        ylabel({'Remove first tau trials','% participants'})
    else
        ylabel('% participants')
    end


end
set(gcf, 'Color', [1 1 1]);
set(gcf, 'PaperUnits', 'centimeters','Units', 'centimeters')
set(gcf,'Position',[0 1 11.6 7])