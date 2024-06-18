function fitPigeon_NLcollapsingbound(dataTable,block_names_publish,num,opt)

arguments
    dataTable    
    block_names_publish
    num = 3
    opt = 'expdecay'
end

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);

switch opt
    case 'expdecay'
        x0=[abs(randn(1).*std(dataTable.bound,'omitnan')+mean(dataTable.bound,'omitnan')),abs(randn(1)),abs(randn(1).*randi(50,1)),randi(50,1)];
        lb = [0,0,0,1];
        ub = [0.75,3,1,51];
        expfcn = @(b,x) b(2).*-exp(b(3).*x-b(4))+b(1);

        fcn = expfcn;
        fdat = nan(numSubjects,4,numBlocks); % 4 fits 
    case 'yparabola'
        x0=abs([randn(1).*std(dataTable.bound,'omitnan')+mean(dataTable.bound,'omitnan'),randn(1),50]);
        lb = [0,0,5];
        ub = [0.75,3,51];
        sqfcn = @(b,x) b(2).*sqrt(-x+b(3))+b(1);
        
        fcn = sqfcn;
        fdat = nan(numSubjects,3,numBlocks); % 3 fits 
end
options = optimset('MaxFunEvals', 10000);

%% Set up figure

figure
tiledlayout(2,3,'TileSpacing','tight','Padding','compact');
EXAMPLE_SUBJECT =13;

%% 

% Loop through each subject
Lg = dataTable.RT>=0;
for ss = 1:numSubjects
    disp(ss)
    for bb = 1:numBlocks
        bounds =[];
        dts =[];
        % Get subject- and block-specific bounds, rts
        Lsb = Lg & dataTable.subjectIndex==subjects(ss) & dataTable.blockIndex==blocks(bb);
        bounds = vertcat(abs(dataTable.boundHi(Lsb)),abs(dataTable.boundLo(Lsb)));
        dts = vertcat(dataTable.DT(Lsb),dataTable.DT(Lsb)-1);
        
        if ss == EXAMPLE_SUBJECT
            egbounds{bb}=bounds;
            egdts{bb} = dts;
        end

        switch opt
            case 'expdecay'
                x0=[abs(randn(1).*std(dataTable.bound(Lsb),'omitnan')+mean(dataTable.bound(Lsb),'omitnan')),abs(randn(1)),abs(randn(1).*randi(50,1)),randi(50,1)];
                
            case 'yparabola'
                x0=abs([randn(1).*std(dataTable.bound(Lsb),'omitnan')+mean(dataTable.bound(Lsb),'omitnan'),randn(1),50]);
               
        end

        % Fit bounds, RTs
        fdat(ss,:,bb) = fmincon(@(b) sum((bounds-fcn(b,dts)).^2), ...
            x0,[],[],[],[],lb,ub,[], options);

    end
end
%% 

for bb = 1:numBlocks
    nexttile(bb); hold on;
    plot(egdts{bb},egbounds{bb},'k.')
    xs = 1:max(egdts{bb});
    plot(xs,fcn(fdat(EXAMPLE_SUBJECT,:,bb),xs),'r')
%     xlim([0 51])
    ylim([0 0.75])
    xlabel('DTs (steps)')
    ylabel('bounds')
    title(block_names_publish(bb))
    
    nexttile(bb+numBlocks)
    switch opt
        case 'expdecay'
            plot(fdat(:,end,bb)+abs((1+fdat(:,end,bb))./fdat(:,3,bb)),fdat(:,1,bb),'ko')
            ylabel('fix start')
            xlabel('decline clear (-exp(1))')
            ylim([0 0.75])
            xlim([0 500])
        case 'yparabola'
            plot(fdat(:,end,bb),fdat(:,1,bb),'ko')
            ylabel('base-y')
            xlabel('base-x')
            xlim([0 51])
            ylim([0 0.75])
    end
    

end
