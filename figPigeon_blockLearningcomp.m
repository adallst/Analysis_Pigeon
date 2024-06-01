
function figPigeon_blockLearningcomp(fdat, fsim, block_names_publish,num)
% function figPigeon_blockLearning(dataTable, num)
%
% Figure: basic performance
%
% top row: RT vs accuracy (% correct), per subject
% bottom row: coins/step vs median bound, per subject

arguments
    fdat    
    fsim
    block_names_publish
    num = 3
end

%% Set up figure
wid     = 17.6; % total width
cols    = {0.65};
hts     = 3.5;

gr = 0.9.*ones(1,3);
wt = 0.99.*ones(1,3);

figure
tiledlayout(3,1)
% [axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 2.5, [], 'Pigeons', true);
% set(axs,'Units','normalized');
for bb = 1:size(fdat,3)
% Right: Scale vs tau
%     axes(axs(bb)); cla reset; hold on;    plot([0 500], [0 0], 'k:')
    nexttile; hold on;
    plot(fdat(:,3,bb), fsim(:,3,bb), 'ko', 'MarkerFaceColor', wt);
%     Lp = fdat(:,4,bb)<0.01;
%     plot(fdat(Lp,3,bb), -fdat(Lp,2,bb), 'ko', 'MarkerFaceColor', 'k');
    medd = nanmedian(fdat(:,3,bb,1));
    meds = nanmedian(fdat(:,3,bb,1));
%     iqrv = iqr(fdat(:,3,bb,1));
%     plot(medd.*[1 1], [-1 1], 'r-')
    yline(nanmedian(medd),'r-')
    yline(nanmedian(meds),'r-')
    plot(0:100,0:100,'k:')
    title(sprintf('med data=%.2f, med sims=%.2f', medd, meds))
    axis([0 200 0 200])
    axis square
    if bb == 3
        xlabel('Data Tau')
        ylabel('Sim Tau')
    end
    
end