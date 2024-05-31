function figPigeon_accuracyVsRT(dataTable, block_names_publish, num, generativeMean, generativeSTD)
% function figPigeon_accuracyVsRT(dataTable, num, generativeMean, generativeSTD)
%
% Figure: accuracy vs RT. duh.
%

arguments
    dataTable  
    block_names_publish
    num = 4
    generativeMean = 0.05
    generativeSTD = 0.15
end

%% Set up figure
wid     = 17.6; % total width
cols    = {3};
hts     = 3.5;
[axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);

RTBins = 1:10;
numRTBins = length(RTBins);
rtData = nan(numSubjects,numRTBins,numBlocks,2); % pct correct, num trials

Lg = dataTable.RT >= 0 & dataTable.correct >= 0;
Lc = dataTable.correct==1;
for bb = 1:numBlocks
    for ss = 1:numSubjects
        for rr = 1:numRTBins
            Lsb = Lg & dataTable.RT==rr & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
            if sum(Lsb)>2
                rtData(ss,rr,bb,1) = sum(Lsb&Lc)./sum(Lsb).*100;
                rtData(ss,rr,bb,2) = sum(Lsb);
            end
        end
    end
end

eData = nan(numRTBins,1);
eData(1) = 50;
for rr = 1:numRTBins-1
    eData(rr+1) = (1-normcdf(0,generativeMean.*rr, sqrt(generativeSTD.^2.*rr))).*100;
end


%% Plotz
% Per block
%gry = 0.9.*ones(1,3);
%wht = 0.99.*ones(1,3);
for bb = 1:numBlocks

    % Set axes
    axes(axs(bb)); cla reset; hold on;
    plot([0 11], [50 50], 'k:');
    
    % Expected RT vs pct correct
    plot(RTBins, eData, 'r-', 'LineWidth', 2);

    % boxplot of RT vs pct correct
    boxplot(rtData(:,:,bb,1), 'labels', RTBins);
%     axis([40 100 0 18])
    if bb == 1
        ylabel('Pct correct')
        xlabel('RT (steps)')
    end
    title(block_names_publish(bb))
end
