function lateslopes = figPigeon_boundByTime4(posall,block_names_publish, ndtbest, params,num)
% function figPigeon_boundByTime(dataTable, num)
%
% Figure: 

arguments
    posall
    block_names_publish
    ndtbest
    params
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
m = params(1);
b = params(2);
mu = params(3);
ml = params(4);



% Collect data per subject/block
subjects = 1:60;
numSubjects = length(subjects);
blocks = 1:3;
numBlocks = length(blocks);

blockStruct = struct;
for bb = 1:numBlocks
    for ss = 1:numSubjects
        bounds(ss,:) = [median(abs(posall{bb,ss}),'omitnan') ...
            std(abs(posall{bb,ss}),'omitnan')];
    end
    blockStruct(bb).blockArgs.blockIndex = bb;
    blockStruct(bb).simArgs.boundMean = bounds(:,1);
    blockStruct(bb).simArgs.boundSTD = 0;%bounds(:,2).*0.01; %most variance comes from slope decrease
    blockStruct(bb).simArgs.boundSlope = m;
    blockStruct(bb).simArgs.boundIntercept = b;
    blockStruct(bb).simArgs.NDTMin = 0;
    blockStruct(bb).simArgs.NDTMax = 0;
    blockStruct(bb).simArgs.lapseRate = 0.00;
end

simDataTablemod = getPigeon_simulatedDataTable2(blockStruct);
dataTable = simDataTablemod;

delays = 0:3;
numDelays = length(delays);
RTs = 2:14;
numRTs = length(RTs);
DTs = 1:14;
numDTs = length(DTs);
bData = nan(numSubjects,numDTs,numBlocks,numDelays,3); % last is mean/sem/n bound
lateslopes = nan(1,4);
earlyslopes = nan(1,4);

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
for dd = 1%:numDelays
    for bb = 2%1:numBlocks

        use =~isnan(bData(:,3:end,bb,1,1));
        thisbdata = bData(:,3:end,bb,1,1);
%         thisRTs = repmat(RTs(4:end),numSubjects,1);
        thisDTs = repmat(DTs(3:end),numSubjects,1);
        mdl2 = fitlm(thisDTs(use),thisbdata(use));
        lateslopes(1,:) = [mdl2.Coefficients.Estimate(2) mdl2.Coefficients.Estimate(1)...
            mdl2.Coefficients.Estimate(2)+mdl2.Coefficients.SE(2) ...
            mdl2.Coefficients.Estimate(2)-mdl2.Coefficients.SE(2)];
       
        % do stats per bin
%         for rr = 1:numDTs
%             Lr = isfinite(bData(:,rr,bb,dd,1));
%             pvals(rr) = signtest(bData(Lr,rr,bb,dd,1),1);
%         end
% 
%         Lp = pvals<0.01./numDTs;          

    end
end
