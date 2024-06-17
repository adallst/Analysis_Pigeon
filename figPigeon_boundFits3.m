function figPigeon_boundFits3(dataTable, block_names_publish, num,name)
% function figPigeon_boundFits(dataTable, num)
%
% Figure: fit bound to linear model
%
% Top panels: Example subject
% Bottom panels: Summary bound fits

arguments
    dataTable
    block_names_publish
    num = 1
    name = "data"
end

%% Set up figure
%
EXAMPLE_SUBJECT = 5;% 13; % 12, 13, 47, 60, 34
wid = 17.6; % total width
hts = 3.5;
cols = {3, 3, 3};
% [axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
% set(axs,'Units','normalized');
% figure
% tiledlayout(3,1)

wt = 0.99.*ones(3,1);
lgr = 0.75.*ones(3,1);

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);
linearFits = nan(numSubjects, 2, numBlocks, 2); % last is fits, slope/intercept

% Set up patternsearch options
patternsearch_opts = optimoptions(@patternsearch,    ...
    'MaxIter',     30000);

% Loop through each subject
Lgood = dataTable.trialNumber>10 & dataTable.RT>3 & dataTable.RT<20;
f0 = [1.3 1.5]; % fit initial values
for ss = 1:numSubjects
    disp(ss)
    for bb = 1:numBlocks

%         % Get subject- and block-specific data
        Lsb = Lgood & dataTable.blockIndex==blocks(bb) & dataTable.subjectIndex==subjects(ss);
        sbdat = cat(2, ...
            abs(dataTable.boundLo(Lsb)), ...
            abs(dataTable.bound(Lsb)), ...
            dataTable.DT(Lsb)-1);
%         medbnd(bb,ss) = median(abs(dataTable.bound(Lsb)),'omitnan');
%         medRT(bb,ss) = median(dataTable.RT(Lsb),'omitnan');
        medDT(bb,ss) = median(dataTable.DT(Lsb),'omitnan');

        % Objective function, wrapper to call boundErr with data
        boundErrFcn = @(fits) boundErrPts(fits, sbdat);

        % Get the fits
        linearFits(ss,:,bb,1) = patternsearch(boundErrFcn, f0, ...
            [],[],[],[],[0.1 0.1],[1.9 1.9],[],patternsearch_opts);

        % compute slope,intercept
        [~,bound,Lint] = boundErrPts(linearFits(ss,:,bb,1), sbdat);
        m = (diff(bound([4 2]))./diff(bound([3 1])));
        b = bound(2)-bound(1)*m;
        linearFits(ss,:,bb,2) = [m b];
        
%         mdl = fitlm(dataTable.DT(Lsb),abs(dataTable.bound(Lsb)));
%         linearFits(ss,:,bb,2) = [mdl.Coefficients.Estimate(2) mdl.Coefficients.Estimate(1)];
%         m=mdl.Coefficients.Estimate(2);
%         b=mdl.Coefficients.Estimate(1);
    

        % Plotz
        if ss==EXAMPLE_SUBJECT &bb==2
%             axes(axs(bb));cla reset; hold on;
            nexttile(num)
            hold on
%             for xx = 1:size(sbdat,1)
%                 h=plot(dataTable.DT(Lsb),abs(dataTable.bound(Lsb)), 'k.');
%                 if Lint(xx)
%                     set(h, 'color', lgr)
%                 end
%             end
%             plot(1:20,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2).*[1:20], 'r-', 'LineWidth', 2)
            for xx = 1:size(sbdat,1)
                h=plot([sbdat(xx,3) sbdat(xx,3)+1], sbdat(xx,1:2), 'k-');
                if Lint(xx)
                    set(h, 'color', lgr)
                end
            end
            plot(bound([1 3]), bound([2 4]), 'r-', 'LineWidth', 2)
            axis([0 20 0 0.8])
            title({name,sprintf('slope=%.2f, y-int=%.2f', linearFits(ss,1,bb,2),linearFits(ss,2,bb,2))})
%             if bb==1
                xlabel('DT (steps)')
                ylabel('Bound magnitude (a.u.)')
%             end
        end
    end
end
    
    %disp(ss)
    %r = input('next')


% Plot summary
    for bb = 2%1:3
    % plot m,b from linear fits, color coded by better model (filled is
    % rt is inv gaussian)  axes(axs(numBlocks+bb))
        nexttile(num+3); cla reset; hold on;    plot([0 0], [-1 1], 'k:');
        plot([-1 1], [0 0], 'k:');
        plot(linearFits(:,2,bb,2), linearFits(:,1,bb,2), 'ko', 'MarkerFaceColor', wt);
        axis([-0.2 0.8 -0.04 0.15])
%         if bb == 1
            xlabel('y-intercept')
            ylabel('Slope')
%         end

        nexttile(num+6); cla reset; hold on;    plot([0 0], [-1 1], 'k:');
    %     plot([-1 1], [0 0], 'k:');axes(axs(numBlocks.*2+bb))
        xline(0,'k:')
        yline(0,'k:')
        plot(medDT(bb,:)', linearFits(:,1,bb,2), 'ko', 'MarkerFaceColor', wt);%medbnd(bb,:)
        nouse = isnan(medDT(bb,:));
        lsline
        [R,P] = corr(medDT(bb,~nouse)', linearFits(~nouse,1,bb,2),'type', 'Pearson');
%         title({block_names_publish(bb),sprintf('rho=%.2f, p=%.3f', R, P)})
        title(sprintf('rho=%.2f, p=%.3f', R, P))
        axis([0 15 -0.04 0.15])
%         if bb == 1
            xlabel('Average DT') %xlabel('avg bound')
            ylabel('Slope')
%         end
    end
% set(gcf, 'Color', [1 1 1]);
% set(gcf, 'PaperUnits', 'centimeters','Units', 'centimeters')
% set(gcf,'Position',[0 1 8.5 17.6])
end
