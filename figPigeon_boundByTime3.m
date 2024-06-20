function lateslopes = figPigeon_boundByTime3(dataTable, block_names_publish, ndtbest, num)
% function figPigeon_boundByTime(dataTable, num)
%
% Figure: 

arguments
    dataTable
    block_names_publish
    ndtbest
    num = 5
end

%% Set up figure
%
% wid     = 17.6; % total width
% cols    = {3,3,3,3};
% hts     = 3.5;
% % [axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
% % set(axs,'Units','normalized');
% figure
% hold on

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);

delays = 0:3;
numDelays = length(delays);
RTs = 2:14;
numRTs = length(RTs);
DTs = 1:14;
numDTs = length(DTs);
bData = nan(numSubjects,numDTs,numBlocks,numDelays,3); % last is mean/sem/n bound
lateslopes = nan(numBlocks,4);
earlyslopes = nan(numBlocks,4);

% Lg = dataTable.RT >= 0 & dataTable.correct >= 0;
Lg = dataTable.DT >= 0 & dataTable.correct >= 0;
for bb = 1:numBlocks
    for ss = 1:numSubjects
        for rr = 1:numDTs
%             Lsb = Lg & dataTable.RT==RTs(rr) & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
            Lsb = Lg & (dataTable.RT-ndtbest(ss,bb))==DTs(rr) & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
            if sum(Lsb)>4
                steps = cell2mat(dataTable.steps(Lsb));
%                 for dd = 1:sum(delays<RTs(rr)-1)
                    bounds = abs(steps(:, max(end-ndtbest(ss,bb),1)));
                    bData(ss,rr,bb,1,:) = [mean(bounds), sem(bounds), length(bounds)];
%                 end
            end
        end
        %normalize to 1?
%         for dd = 1:numDelays
            bData(ss,:,bb,1,1) = bData(ss,:,bb,1,1)./ ...
                (sum(bData(ss,:,bb,1,1).*bData(ss,:,bb,1,3),'omitnan')./ ...
                sum(bData(ss,:,bb,1,3),'omitnan'));
            
%         end
    end
end


%% Plotz
% Per block
%gry = 0.9.*ones(1,3);
%wht = 0.99.*ones(1,3);
pvals = nan(numDTs,1);
hold on
for dd = 1%:numDelays
    for bb = 2%1:numBlocks
%          use =~isnan(bData(:,1:5,bb,1,1));
%         thisbdata = bData(:,1:5,bb,1,1);
% %         thisRTs = repmat(RTs(4:end),numSubjects,1);
%         thisDTs = repmat(DTs(1:5),numSubjects,1);
%         mdl1 = fitlm(thisDTs(use),thisbdata(use));
%         earlyslopes(bb,:) = [mdl1.Coefficients.Estimate(2) mdl1.Coefficients.Estimate(1)...
%             mdl1.Coefficients.Estimate(2)+mdl1.Coefficients.SE(2) ...
%             mdl1.Coefficients.Estimate(2)-mdl1.Coefficients.SE(2)];
%         % Set axes
% %         axes(axs((dd-1)*numBlocks+bb)); cla reset; hold on;
%         plot([0 11], [0.5 0.5], 'k:');
%         plot([0 11], [1 1], 'k:');
%         yline(0,'k:')
% 
%         % boxplot of RT vs bound bin
% %         boxplot(bData(:,:,bb,1,1), 'labels', RTs);
%         ad_boxplot(bData(:,:,bb,1,1),'skyline',0)
%         xticks(1:size(bData(:,:,bb,1,1),2))
%         xticklabels(DTs(1:size(bData(:,:,bb,1,1),2)))
%         ylim([-0.05 2]);
% 
%         % do stats per bin
%         for rr = 1:numDTs
%             Lr = isfinite(bData(:,rr,bb,dd,1));
%             pvals(rr) = signtest(bData(Lr,rr,bb,dd,1),1);
%         end
% 
%         Lp = pvals<0.01./numDTs;
% %         plot(find(Lp), 2.*ones(sum(Lp),1), 'mo', 'MarkerFaceColor', 'm');
%         
%         plot(DTs(1:5),earlyslopes(bb,1).*DTs(1:5)+earlyslopes(bb,2),'r')
% 
%         if bb == 1
%             xlabel('DT (steps)')
%             ylabel(['NDT = ' num2str(delays(dd))])
%         elseif bb>1
%             ylabel('Bounds rel. mean')
%             xlabel('DT (steps)')
%         end
%         if dd == 1
% %             title({num,block_names_publish(bb)})
%             title(num)
%         end
%         mdl = fitlm(RTs(4:end),nanmedian(bData(:,4:end,bb,1,1)));
%         lateslopes(bb,:) = [mdl.Coefficients.Estimate(2) mdl.Coefficients.Estimate(1)...
%             mdl.Coefficients.Estimate(2)+mdl.Coefficients.SE(2) ...
%             mdl.Coefficients.Estimate(2)-mdl.Coefficients.SE(2)];
%         nexttile; hold on
        use =~isnan(bData(:,3:end,bb,1,1));
        thisbdata = bData(:,3:end,bb,1,1);
%         thisRTs = repmat(RTs(4:end),numSubjects,1);
        thisDTs = repmat(DTs(3:end),numSubjects,1);
        mdl2 = fitlm(thisDTs(use),thisbdata(use));
        lateslopes(bb,:) = [mdl2.Coefficients.Estimate(2) mdl2.Coefficients.Estimate(1)...
            mdl2.Coefficients.Estimate(2)+mdl2.Coefficients.SE(2) ...
            mdl2.Coefficients.Estimate(2)-mdl2.Coefficients.SE(2)];
        % Set axes
%         axes(axs((dd-1)*numBlocks+bb)); cla reset; hold on;
        plot([0 11], [0.5 0.5], 'k:');
        plot([0 11], [1 1], 'k:');
        yline(0,'k:')

        % boxplot of RT vs bound bin
%         boxplot(bData(:,:,bb,1,1), 'labels', RTs);
        ad_boxplot(bData(:,:,bb,1,1),'skyline',0)
        xticks(0:4:size(bData(:,:,bb,1,1),2))
        xticklabels(0:4:size(bData(:,:,bb,1,1),2))
        ylim([-0.05 2]);

        % do stats per bin
        for rr = 1:numDTs
            Lr = isfinite(bData(:,rr,bb,dd,1));
            pvals(rr) = signtest(bData(Lr,rr,bb,dd,1),1);
        end

        Lp = pvals<0.01./numDTs;
%         plot(find(Lp), 2.*ones(sum(Lp),1), 'mo', 'MarkerFaceColor', 'm');
        
        plot(DTs(3:end),lateslopes(bb,1).*DTs(3:end)+lateslopes(bb,2),'r')

        if bb == 1
            xlabel('DT (steps)')
            ylabel(['NDT = ' num2str(delays(dd))])
        elseif bb>1
            ylabel('Bounds rel. mean')
            xlabel('DT (steps)')
        end
        if dd == 1
%             title({num,block_names_publish(bb)})
            title(num)
        end
        
%         nexttile; hold on
%         plot(repmat(DTs(1:size(bData(:,:,bb,1,1),2)),size(bData(:,:,bb,1,1),1),1),bData(:,:,bb,1,1),'k.')
%         plot(DTs(1:size(bData(:,:,bb,1,1),2)),median(bData(:,:,bb,1,1),'omitnan'),'k')
%         ylim([-0.05 4]);
%         ylabel('Bounds rel. mean')
%         xlabel('DT (steps)')
%         title(num)
    end
end
set(gcf, 'Color', [1 1 1]);
set(gcf, 'PaperUnits', 'centimeters','Units', 'centimeters')
set(gcf,'Position',[0 2 8.5 8.5])
