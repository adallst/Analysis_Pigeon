function fdat = figPigeon_blockLearning2(dataTable, block_names_publish, num)
% function figPigeon_blockLearning(dataTable, num)
%
% Figure: basic performance
%
% top row: RT vs accuracy (% correct), per subject
% bottom row: coins/step vs median bound, per subject

arguments
    dataTable    
    block_names_publish
    num = 3
end

%% Set up figure
wid     = 17.6./3; % total width
cols    = {[0.65 0.25]};
hts     = [3.5];
% [axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 2.5, [], 'Pigeons', true);
% set(axs,'Units','normalized');
figure
tiledlayout(3,4,'TileSpacing','tight','Padding','compact');
EXAMPLE_SUBJECT =12;
% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);

% Count median, max trials for data matrix
maxTrials = 0;
Lg = dataTable.RT>=0;
for ss = 1:numSubjects
    for bb = 1:numBlocks
        Lsb = Lg & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
        if sum(Lsb) > maxTrials
            maxTrials = sum(Lsb);
        end
    end
end

% For exponential fits
fcn = @(b,x) b(1)+b(2).*exp(-x./b(3));
options = optimset('MaxFunEvals', 10000);
fdat = nan(numSubjects,4,numBlocks); % 2: fits + pvalue; 4: bound/RT

% Loop through each subject
zBounds = nan(numSubjects,maxTrials,numBlocks);
for ss = 1:numSubjects
    disp(ss)
    for bb = 1:numBlocks

        % Get subject- and block-specific bounds, rts
        Lsb = Lg & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
        bounds = abs(dataTable.bound(Lsb));
        % rts = dataTable.RT(Lsb);
        xs = (1:sum(Lsb))';

        % compare start/end bounds
        fdat(ss,4,bb) = ranksum(bounds(1:10), bounds(end-9:end));

        % Fit bounds, RTs
        fdat(ss,1:3,bb) = fmincon(@(b) sum((bounds-fcn(b,xs)).^2), ...
            [mean(bounds(1:10)) mean(bounds(end-10:end))./mean(bounds(1:10)) 20],[],[],[],[],[0 -0.75 1],[0.75 0.75 100],[], options);
        %fdat(ss,:,bb,2) = fmincon(@(b) sum((rts-fcn(b,xs)).^2), ...
        %    [0 0.75 20],[],[],[],[],[0 -20 1],[20 20 100],[], options);
        %
        %         subplot(2,1,1); cla reset; hold on;
        %         plot(xs, bounds, 'ko');
        %         plot(xs, fcn(fdat(ss,:,bb,1),xs), 'r-')
        %         title(sprintf('Lower=%.2f, Upper=%.2f, Tau=%.2f', ...
        %             fdat(ss,1,bb,1), fdat(ss,2,bb,1), fdat(ss,3,bb,1)))
        %
        %         subplot(2,1,2); cla reset; hold on;
        %         plot(xs, rts, 'ko');
        %         plot(xs, fcn(fdat(ss,:,bb,2),xs), 'r-')
        %         title(sprintf('Lower=%.2f, Upper=%.2f, Tau=%.2f', ...
        %             fdat(ss,1,bb,2), fdat(ss,2,bb,2), fdat(ss,3,bb,2)))
        %
        %         r = input('next')
        
        if ss==EXAMPLE_SUBJECT
            nexttile(3+4.*(bb-1)); hold on;
            plot(xs, bounds, 'k.');
            plot(xs, fcn(fdat(ss,:,bb),xs), 'r-',"LineWidth",2)
%             title(sprintf('Lower=%.2f, Upper=%.2f, Tau=%.2f', ...
%                 fdat(ss,1,bb,1), fdat(ss,2,bb,1), fdat(ss,3,bb,1)))
            title("Example Fit")
            xlabel("Trials")
            ylabel("Bound (a.u.)")
        end
        

        % z-score bound from full block
        zBounds(ss,1:length(bounds),bb) = zscore(bounds);
    end
end
%% Regression to mean
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
subjectCutoff = 0.5; % exclude trials with less than this frac of subj data

for bb = 1:numBlocks

    % Left: zbound median, IQR
    nexttile((1+4.*(bb-1)),[1,2]); hold on;
    subjectFraction = sum(isfinite(zBounds(:,:,bb)))./numSubjects;
    trialCutoff = find(subjectFraction<subjectCutoff,1);
    xax = 1:trialCutoff;
    pcts = prctile(zBounds(:,1:trialCutoff,bb), [25 50 75], 1);
    xp = [xax flip(xax)];
    yp = [pcts(1,:) flip(pcts(3,:))];
    Lp = isfinite(yp);
    h=patch(xp(Lp), yp(Lp), gr);
    set(h, 'LineStyle', 'none')
    plot(xax, pcts(2,:), 'k-', 'LineWidth', 2);
    plot([1 120], [0 0], 'k:');
    axis([1 120 -1.5 1.5])
%     if bb == 3
        xlabel('Trial number')
%         ylabel('Bound (z-score)')
%     end
    ylabel({block_names_publish(bb),'Bound (z-score)'})
    
    
    


    % Right: Scale vs tau
    nexttile(4+4.*(bb-1)); hold on;
    plot([0 500], [0 0], 'k:')
    plot(fdat(:,3,bb), -fdat(:,2,bb), 'ko', 'MarkerFaceColor', wt);
    Lp = fdat(:,4,bb)<0.01;
    plot(fdat(Lp,3,bb), -fdat(Lp,2,bb), 'ko', 'MarkerFaceColor', 'k');
    medv = median(fdat(:,3,bb,1));
    iqrv = iqr(fdat(:,3,bb,1));
    meds = median(-fdat(Lp,2,bb));
    plot(medv.*[1 1], [-1 1], 'r-')
    yline(median(-fdat(:,2,bb,1)),'r-')
%     title(sprintf('med=%.2f, iqr=%.2f', medv, iqrv))
%     title(sprintf('med tau=%.2f \n med scale=%.2f', medv, meds),'FontSize',8)
    axis([0 100 -1 1])
%     if bb == 3
        xlabel('Tau')
        ylabel('Scale')
%     end
    
%     nexttile(9+bb);hold on
%     ct =1;
%     bins = -1.5:0.125:0.1;
% %     plot(rData(:,1,bb,ct), rData(:,2,bb,ct), 'ko', 'MarkerFaceColor', wt);
%     Lsig = [];
%     Lsig = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
% %     plot(rData(Lsig,1,bb,ct), rData(Lsig,2,bb,ct), 'ko', 'MarkerFaceColor', 'k');
% %     plot([0 0.75], [-1 -1], 'k:');
% %     axis([0 0.75 -2 1])
%     h1 = histcounts(rData(Lsig,2,bb,ct),'BinEdges',bins);
%     h2 = histcounts(rData(~Lsig,2,bb,ct),'BinEdges',bins);
%     bar(bins(2:end),[h1;h2]'./numSubjects,'stacked')
%     newcolors = [0 0 0; 1 1 1];
%     colororder(newcolors)
% %     if bb == 3
% %         xlabel('Intercept')
% %         ylabel('Slope')
%         xlabel('Slope')
% %     end
%     if bb ==1
%         ylabel({['Remove first ' num2str(rmhead(ct)) ' trials'],'% participants'})
%     else
%         ylabel('% participants')
%     end
%     title(block_names_publish(bb))
%     
%     nexttile(12+bb); hold on;
%     ct =2;
% %     plot(rData(:,1,bb,ct), rData(:,2,bb,ct), 'ko', 'MarkerFaceColor', wt);
%     Lsig = [];
%     Lsig = rData(:,3,bb,ct)<0.005 | rData(:,3,bb,ct)>0.995;
% %     plot(rData(Lsig,1,bb,ct), rData(Lsig,2,bb,ct), 'ko', 'MarkerFaceColor', 'k');
% %     plot([0 0.75], [-1 -1], 'k:');
% %     axis([0 0.75 -2 1])
%     h1 = histcounts(rData(Lsig,2,bb,ct),'BinEdges',bins);
%     h2 = histcounts(rData(~Lsig,2,bb,ct),'BinEdges',bins);
%     bar(bins(2:end),[h1;h2]'./numSubjects,'stacked')
%     newcolors = [0 0 0; 1 1 1];
%     colororder(newcolors)
% %     if bb == 3
% %         xlabel('Intercept')
% %         ylabel('Slope')
%         xlabel('Slope')
% %     end
%     if bb ==1
%         ylabel({'Remove first tau trials','% participants'})
%     else
%         ylabel('% participants')
%     end
    
%     nexttile(9+bb);hold on
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
end
set(gcf, 'Color', [1 1 1]);
set(gcf, 'PaperUnits', 'centimeters','Units', 'centimeters')
set(gcf,'Position',[0 1 15.6 11.6])
