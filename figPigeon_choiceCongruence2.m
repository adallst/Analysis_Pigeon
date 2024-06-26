function ndtbest = figPigeon_choiceCongruence2(dataTable, block_names_publish, num)
% function figPigeon_choiceCongruence(dataTable, num)
%
% Figure:
%

arguments
    dataTable
    block_names_publish
    num = data
end

%% Set up figure
wid     = 17.6; % total width
cols    = {3,3,3,3};
hts     = 3.5;
% [axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
% set(axs,'Units','normalized');
figure
tiledlayout(6,2,'TileSpacing','compact')

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
ndtbest = nan(numSubjects,numBlocks);
medRT = nan(numSubjects,numBlocks);

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
        if sum(sum(isnan(ccData(ss,4:end,bb,:)))) ~= numel(ccData(ss,4:end,bb,:))
        ndtbest(ss,bb)=find(nanmean(ccData(ss,4:end,bb,:))==max(nanmean(ccData(ss,4:end,bb,:))),1)-1;%starts at 1, change to 0
        medRT(ss,bb) = median(dataTable.RT(Lg & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb)));
        end
    end
end

ndtmedRT = nan(50,numDelays);

%% 
% figPigeon_ndt(dataTable, simDataTable, 1.1)


%% Plotz
% Per block
gry = 0.9.*ones(1,3);
%wht = 0.99.*ones(1,3);
% firstcol = [1:4:13];
for dd = 1:numDelays
    for bb = 2%1:numBlocks

        % Set axes
        nexttile([2,1]); cla reset; hold on;
        plot([0 11], [0.5 0.5], 'k:');
        plot([0 11], [1 1], 'k:');

        % boxplot of RT vs congruence
%         boxplot(ccData(:,:,bb,dd), 'labels', RTs,'symbol','.','Colors',[0.2 0.2 0.2]);
        ad_boxplot(ccData(:,:,bb,dd),'rowhouse',0)
        xticks(1:min(factor(size(ccData(:,:,bb,dd),2))):size(ccData(:,:,bb,dd),2))
        xticklabels(RTs(1:min(factor(size(ccData(:,:,bb,dd),2))):size(ccData(:,:,bb,dd),2)))
        ylim([0 1.5])

        %     axis([40 100 0 18])
        if dd == 1|dd==3
            xlabel('RT (steps)')
%             ylabel(['NDT = ' num2str(delays(dd))])
            ylabel('Congruence')
        elseif dd~=3|dd~=1
%             ylabel('Congruence')
            xlabel('RT (steps)')
        end
%         if dd == 1
            title(['NDT = ' num2str(delays(dd))])
%         end
        
        ndtmedRT(1:length(medRT(ndtbest(:,bb)==delays(dd),bb)),dd) = medRT(ndtbest(:,bb)==delays(dd),2);
    end
end
nexttile(9,[2 2])
hold on
% boxplot(ndtmedRT,'labels',delays,'Colors','k','notch', 'on')
% ad_boxplot(ndtmedRT,'gray')
% s = swarmchart(repmat(delays+1,50,1),ndtmedRT,4,'k',XJitterWidth = 0.5);
% plot(delays+1,median(ndtmedRT,'omitnan'),'r_',"MarkerSize",12)
% xticks(1:length(delays))
% xticklabels(arrayfun(@num2str,delays,'UniformOutput',false))
% plot(repmat(delays+1,length(ndtmedRT),1)+randn(size(ndtmedRT)).*0.075,ndtmedRT,'.k')
histogram(ndtbest(:,2),"FaceColor",gry)
xlabel('NDT')
ylabel('# participants')
% ylabel('Median RT (steps)')
% sgtitle(num)
set(gcf, 'Color', [1 1 1]);
set(gcf, 'PaperUnits', 'centimeters','Units', 'centimeters')
set(gcf,'Position',[0 0 6 10])
set(gcf,'PaperPosition',[0 0 6 10])