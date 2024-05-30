function figPigeon_choiceCongruence(dataTable, num)
% function figPigeon_choiceCongruence(dataTable, num)
%
% Figure:
%

arguments
    dataTable
    num = 4
end

%% Set up figure
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
RTs = 2:10;
numRTs = length(RTs);
ccData = nan(numSubjects,numRTs,numBlocks,numDelays);

Lg = dataTable.RT >= 0 & dataTable.correct >= 0;
for bb = 1:numBlocks
    for ss = 1:numSubjects
        for rr = 1:numRTs
            Lsb = Lg & dataTable.RT==RTs(rr) & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
            if any(Lsb)
                steps = cell2mat(dataTable.steps(Lsb));
                choices = dataTable.choice(Lsb)-0.5;
                for dd = 1:sum(delays<RTs(rr))
                    ccData(ss,rr,bb,dd) = sum(sign(steps(:,end-delays(dd))) == sign(choices))./length(choices);
                end
            end
        end
    end
end


%% Plotz
% Per block
%gry = 0.9.*ones(1,3);
%wht = 0.99.*ones(1,3);
for dd = 1:numDelays
    for bb = 1:numBlocks

        % Set axes
        axes(axs((dd-1)*numBlocks+bb)); cla reset; hold on;
        plot([0 11], [0.5 0.5], 'k:');
        plot([0 11], [1 1], 'k:');

        % boxplot of RT vs pct correct
        boxplot(ccData(:,:,bb,dd), 'labels', RTs);
        %     axis([40 100 0 18])
        %     if bb == 1
        %         xlabel('Pct correct')
        %         ylabel('RT (steps)')
        %     end
    end
end