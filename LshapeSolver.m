function out = LshapeSolver(stageProbs,rpm,opts)
% Lshaped method for Two-Stage Stochastic LP
% Inputs:
% stageProbs{t}: the main program at stage t
% rpm: uncertainty at stage t
% opts: algorithm options
% - opts.x1: initial iterate at root stage
% - opts.verbose: 1 for printing the iteration information
% - opts.savedir/logdir: iterate log save path
% - opts.samplepath: the provided sampling paths if any
% - opts.maxIt: the maximal iteration

% Author: Zhiyuan Zhang
% Date: May 25, 2025


if ~isfield(opts,"maxIt")
    fprintf("[L-shaped] Set the maximal iteration to 100.\n");
    maxIt = 100;
else
    maxIt = opts.maxIt;
end


if ~isfield(opts,"tolerance")
    opts.tolerance = 1e-6;
end

if ~isfield(opts,"x1")
    fprintf("[L-shaped] No initial point is provided.\n");
    x1 = [];
else
    x1 = opts.x1;
end

% prepare for the single cut variable [eta]
for t = 1 : length(stageProbs)-1
    nt = size(stageProbs{t}.D,2);
    stageProbs{t}.D(:,end+1) = 0;
    stageProbs{t}.varnames{end+1} = 'eta';
    stageProbs{t}.beta = sparse(nt,0);
    stageProbs{t}.neg_beta_eta = sparse(nt+1,0);
    stageProbs{t}.alpha = [];
    if isfield(stageProbs{t},'ub')
        stageProbs{t}.ub(end+1) = inf;
    end
    if isfield(stageProbs{t},'lb')
        stageProbs{t}.lb(end+1) = -inf;
    end
end

% loginfo = [];
% loginfo.ocnum = zeros(length(stageProbs),1);
% loginfo.fcnum = zeros(length(stageProbs),1);

out = [];
out.ub = zeros(maxIt,1);
out.lb = zeros(maxIt,1);

obj_exact = zeros(maxIt,1);
for k = 1 : maxIt
    [stageProbs,outk] =  LshapeCore(stageProbs,rpm,opts,k,x1);
    % [stageProbs,out] =  LshapeCoreQP(stageProbs,rpm,opts,k,x1);
    obj_exact(k,1) = outk.lb;
    if outk.stop
        out.iter = k;
        % x_iters = x_iters(:,1:k);
        obj_exact = obj_exact(1:k,1);
        out.x = outk.x;
        break
    end
    x1 = [];
end
out.iter = k;
out.lb = obj_exact;
end