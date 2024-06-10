function [choices, rts, DV, bounds, ndts] = getPigeon_simulatedData(options)
%
% Returns:
%   choices ... 0/1
%   rts     ... values (really indices into DV)
%   DV      ... decision variable (step matrix, rows are trials, cols are
%                   steps)
%   bounds  ... array of bounds per trial
%   ndts    ... array of ndts per trial
%
% Check arguments
arguments
    options.generativeMean      double = 0.05   % generative mean
    options.generativeSTD       double = 0.15   % generative std
    options.numTrials           double = 100000 % Simulated trials per SNR
    options.maxStepsPerTrial    double = 50     % Total possible samples per "trial"
    options.boundMean           double = 0      % bound mean/start value
    options.boundSTD            double = 0      % bound variability (trial-to-trial)
    options.boundSlope          double = 0      % Time-dependent slope
    options.boundMax            double = 0.75;  % max, fixed by task geometry
    options.NDTMin              double = 0;     % non-decision time
    options.NDTmax              double = 0;     % non-decision time
    options.lapseRate           double = 0;     % fraction "guess" trials
end

% Make the DV per trial
DV = cat(2, zeros(options.numTrials, 1), cumsum(normrnd(...
    options.generativeMean, ...
    options.generativeSTD, ...
    options.numTrials, ...
    options.maxStepsPerTrial),2));

% Make matrix of bounds per trial per time
if options.boundSTD > 0
    boundPerTrial = max(0, normrnd(options.boundMean, options.boundSTD, options.numTrials, 1));
else
    boundPerTrial = repmat(options.boundMean, options.numTrials, 1);
end
boundMatrix = repmat(boundPerTrial, 1, options.maxStepsPerTrial);
if options.boundSlope~=0
%     boundMatrix = boundMatrix .* repmat(cat(2,linspace(1,options.boundSlope,10), ...
%         options.boundSlope.*ones(1,options.maxStepsPerTrial-10)), options.numTrials, 1);
    boundMatrix = boundMatrix .* repmat((1:options.boundSlope:(0+abs(options.boundSlope)))+0.5, options.numTrials, 1);
end
boundMatrix(boundMatrix>options.boundMax)=options.boundMax;

% Make array of non-decision times
if options.NDTmax > options.NDTMin
%     ndts = randi(options.NDTmax-options.NDTMin+1,options.numTrials,1)-options.NDTMin+1;
    if options.NDTMin ==0
        ndts = randi((options.NDTmax-options.NDTMin)+1,options.numTrials,1)-1;
    else
        ndts = randi(options.NDTmax-options.NDTMin,options.numTrials,1)+options.NDTMin;
    end
else
    if options.NDTMin ==0
        ndts = options.NDTMin*zeros(options.numTrials,1);
    else
        ndts = options.NDTMin*ones(options.numTrials,1);
    end
end

% To save
choices = nans(options.numTrials, 1);
rts = nans(options.numTrials, 1);
bounds = nans(options.numTrials, 1);

% disp(sprintf('looping through %d trials', options.numTrials))
% Loop through the trials to find first crossings

aDV = abs(DV);
Lcrossed = false(options.numTrials,1);
for ss = 1:options.maxStepsPerTrial
    Lrt = aDV(:,ss)>=boundMatrix(:,ss);
    if any(~Lcrossed & Lrt)
        Lrt = ~Lcrossed & Lrt;
        bounds(Lrt) = boundMatrix(Lrt,ss);
        Lcrossed = Lcrossed | Lrt;
        if ss == 1
            rts(Lrt) = ss+randi(2,sum(Lrt),1)-1;
            choices(Lrt) = randi(2,sum(Lrt),1)-1;
        else
            rts(Lrt) = ss+ndts(Lrt);
            choices(Lrt) = double(DV(Lrt,ss)>0);
        end
    end
end
%LtinyBound = bounds<0.05;
%rts(LtinyBound) = rts(LtinyBound) - ndts(LtinyBound);
%fprintf('median bound in=%.2f, out=%.2f\n', options.boundMean, median(bounds,'omitnan'))

% for tt = 1:options.numTrials-numLapses
%     thisbound = tbound.*bound(tt);
%     rt = find(aDV(tt,:)>=thisbound, 1);
%     if ~isempty(rt)
%         choices(tt) = double(DV(tt,rt)>0);
%         bounds(tt) = thisbound(rt);
%         rts(tt) = rt+ndts(tt);
%     end
% end
rts(~isfinite(rts)) = options.maxStepsPerTrial;
rts(rts>options.maxStepsPerTrial) = options.maxStepsPerTrial;

% Add lapses
if options.lapseRate > 0
    Llapse = false(options.numTrials,1);
    Llapse(randperm(options.numTrials, ceil(options.lapseRate*options.numTrials))) = true;
    rts(Llapse) = randi(15, sum(Llapse), 1);
    choices(Llapse) = randi(2, sum(Llapse), 1)-1;
    bounds(Llapse) = boundMatrix(sub2ind(size(boundMatrix), find(Llapse), rts(Llapse)));
end

% % Loop through the trials to find first crossings
% if options.boundSlope==0
%     % Fixed bound
%     for tt = 1:options.numTrials
%         rt = find(DV(tt,:)>=bound(tt) | DV(tt,:)<=-bound(tt), 1);
%         if ~isempty(rt)
%             rt = min(rt+ndts(tt), options.maxStepsPerTrial);
%             choices(tt) = double(DV(tt,rt)>0);
%             rts(tt) = rt;
%             steps{tt} = DV(tt,1:rt);
%             bounds(tt) = bound(tt);
%         end
%     end
% else
%     % Time-dependent bound
%     %    tbound = cat(2,linspace(1,options.boundSlope,10), linspace(options.boundSlope,0,options.maxStepsPerTrial-10));
%     tbound = cat(2,linspace(1,options.boundSlope,10), options.boundSlope.*ones(1,options.maxStepsPerTrial-10));
%     for tt = 1:options.numTrials
%         thisbound = tbound.*bound(tt);
%         rt = find(DV(tt,:)>=thisbound | DV(tt,:)<=-thisbound, 1);
%         if ~isempty(rt)
%             rt = min(rt+ndts(tt), options.maxStepsPerTrial);
%             choices(tt) = double(DV(tt,rt)>0);
%             rts(tt) = rt;
%             steps{tt} = DV(tt,1:rt);
%             bounds(tt) = thisbound(rt);
%         end
%     end
% end
%
% % guess on lapse trials
% if options.lapseRate > 0
%     lapses = randperm(options.numTrials, round(options.lapseRate*options.numTrials));
%     choices(lapses) = -choices(lapses); %2*randi(2, length(lapses), 1)-3;
% end
