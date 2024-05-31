function figPigeon_boundByTime(dataTable, block_names_publish, num)
% function figPigeon_boundByTime(dataTable, num)
%
% Figure: 

arguments
    dataTable
    block_names_publish
    num = 5
end

%% Set up figure
%
wid     = 17.6; % total width
cols    = {3,3,3,3};
hts     = 3.5;
[axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);

delays = 0:3;
numDelays = length(delays);
RTs = 2:14;
numRTs = length(RTs);
bData = nan(numSubjects,numRTs,numBlocks,numDelays,3); % last is mean/sem/n bound

Lg = dataTable.RT >= 0 & dataTable.correct >= 0;
for bb = 1:numBlocks
    for ss = 1:numSubjects
        for rr = 1:numRTs
            Lsb = Lg & dataTable.RT==RTs(rr) & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
            if sum(Lsb)>4
                steps = cell2mat(dataTable.steps(Lsb));
                for dd = 1:sum(delays<RTs(rr)-1)
                    bounds = abs(mean(steps(:, end-delays(dd)-1:end-delays(dd)),2));
                    bData(ss,rr,bb,dd,:) = [mean(bounds), sem(bounds), length(bounds)];
                end
            end
        end
        %normalize to 1?
        for dd = 1:numDelays
            bData(ss,:,bb,dd,1) = bData(ss,:,bb,dd,1)./ ...
                (sum(bData(ss,:,bb,dd,1).*bData(ss,:,bb,dd,3),'omitnan')./ ...
                sum(bData(ss,:,bb,dd,3),'omitnan'));
        end
    end
end


%% Plotz
% Per block
%gry = 0.9.*ones(1,3);
%wht = 0.99.*ones(1,3);
pvals = nan(numRTs,1);
for dd = 1:numDelays
    for bb = 1:numBlocks

        % Set axes
        axes(axs((dd-1)*numBlocks+bb)); cla reset; hold on;
        plot([0 11], [0.5 0.5], 'k:');
        plot([0 11], [1 1], 'k:');

        % boxplot of RT vs bound bin
        boxplot(bData(:,:,bb,dd,1), 'labels', RTs);
        ylim([0 2]);

        % do stats per bin
        for rr = 1:numRTs
            Lr = isfinite(bData(:,rr,bb,dd,1));
            pvals(rr) = signtest(bData(Lr,rr,bb,dd,1),1);
        end

        Lp = pvals<0.01./numRTs;
        plot(find(Lp), 2.*ones(sum(Lp),1), 'mo', 'MarkerFaceColor', 'm');

        if bb == 1
            xlabel('RT (steps)')
            ylabel(['NDT = ' num2str(delays(dd))])
        elseif bb>1
            ylabel('bounds rel. mean')
            xlabel('RT (steps)')
        end
        if dd == 1
            title(block_names_publish(bb))
        end
    end
end
