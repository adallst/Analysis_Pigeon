function figPigeon_boundTertiles(dataTable, block_names_publish, num)
% function figPigeon_boundTertiles(dataTable, num)
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
cols    = {3,3};
hts     = 3.5;
[axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);

numRTBins = 3;
bData = nan(numSubjects,numRTBins,numBlocks);
bounds = abs(dataTable.bound);

Lg = dataTable.trialNumber > 0 & dataTable.RT >= 0 & dataTable.RT <= 10 & dataTable.correct >= 0;
for bb = 1:numBlocks
    for ss = 1:numSubjects
        Lsb = Lg & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
        rtBins = prctile(dataTable.RT(Lsb), [0 25 75 100]); %linspace(0,100,numRTBins+1));
        for ii = 1:numRTBins
            Lb = Lsb & dataTable.RT > rtBins(ii) & dataTable.RT <= rtBins(ii+1);
            if sum(Lsb)>3
                bData(ss,ii,bb) = median(abs(bounds(Lb)), 'omitnan');
            end
        end
    end
end

%% Plotz
% Per block
%gry = 0.9.*ones(1,3);
%wht = 0.99.*ones(1,3);
for bb = 1:numBlocks
    for ii = 1:2 % early vs mid, late vs mid

        % Set axes
        axes(axs((ii-1)*numBlocks+bb)); cla reset; hold on;
        plot([0 0.75], [0 0.75], 'k:');

        % scatterplot of mean bound in different RT tertiles
        xs = bData(:,2,bb);
        ys = bData(:,(ii-1)*2+1,bb);
        Lg = isfinite(xs) & isfinite(ys);
        plot(xs, ys, 'ko');
        P = signrank(xs(Lg), ys(Lg)); 
%         title(sprintf('p=%.2f', P))
        title({block_names_publish(bb),sprintf('WSR p=%.2f', P)}) 

        axis([0 0.8 0 0.8])
        if bb == 1 && ii == 1
            xlabel('Mean bound mid RT')
            ylabel('Mean bound fast RT')
        elseif bb == 1 && ii == 2
            xlabel('Mean bound mid RT')
            ylabel('Mean bound slow RT')
        end
%         title(block_names_publish(bb))
    end
end
