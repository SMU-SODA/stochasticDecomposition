function out = sdSolver(stageProbs,rpm,opts)
% Stochastic Decomposition for Two-Stage Stochastic LP
% Inputs:
% stageProbs{t}: the main program at stage t
% rpm: uncertainty at stage t
% opts: algorithm options
% - opts.x1: initial iterate at root stage
% - opts.constL: the lower bound of cost-to-go function
% - opts.verbose: 1 for printing the iteration information
% - opts.savedir/logdir: iterate log save path
% - opts.samplepath: the provided sampling paths if any
% - opts.maxIt: the maximal iteration
% - opts.writeflag: 1 for exporting the intermediate programs
% Author: Zhiyuan Zhang
% Date: May 25, 2025


if ~isfield(opts,"maxIt")
    fprintf("[SD] Set the maximal iteration to 100 by default.\n");
    maxIt = 100;
else
    maxIt = opts.maxIt;
end

if ~isfield(opts,"x1")
    fprintf("[SD] No initial point is provided.\n");
    x1 = [];
else
    x1 = opts.x1;
end

if ~isfield(opts,"writeflag")
    fprintf("[SD] Do not export the programs.\n");
    opts.writeflag = 0;
end

if ~isfield(opts,"constL")
    fprintf("[SD] Set the lower bound of the recourse to 0 by default.\n");
    constL = 0;
else
    constL = opts.constL;
end



if ~isfield(opts,"tolerance")
    opts.tolerance = 1e-8;
end

if ~isfield(opts,"savedir")
    opts.savedir = './sdlogs';
end

if ~isfield(opts,"verbose")
    opts.verbose = 0;
end

if ~isfield(opts,"logdir")
    opts.logdir = fullfile(opts.savedir,'default');
end

if ~isfolder(opts.logdir)
    mkdir(opts.logdir)
end

% prepare for the single cut variable [eta]
for t = 1 : length(stageProbs)-1
    nt = size(stageProbs{t}.D,2);
    stageProbs{t}.beta = sparse(nt,0);
    stageProbs{t}.beta_inc = sparse(nt,0);
    stageProbs{t}.alpha = [];
    stageProbs{t}.alpha_inc = [];
end

out = [];
out.ub = zeros(maxIt,1);
out.lb = zeros(maxIt,1);


bundlePI = [];
bundlePI.size = 0;
bundlePI.pi = [];
bundlePI.dualbd = [];
bundlePI.iter = [];

bundleJ = [];
bundleJ.indices = [];
bundleJ.iteridx = [];
bundleJ.alpha = [];
bundleJ.beta = [];
bundleJ.cnt = 0;

OmegaInfo = [];
OmegaInfo.indices = [];
OmegaInfo.cnt = 0;
OmegaInfo.weights = [];

%initial stage
for k = 1 : maxIt
    % sample = [];
    % load samplepath or randomly generate a sample path
    if isfield(opts,"samplepath")
        observ = opts.samplepath(k,:);
    else
        val = rand;
        cumm = 0;
        select_idx = rpm.S;
        for ii = 1 : rpm.S
            cumm = cumm + rpm.probs(ii);
            if cumm > val
                select_idx = ii;
                break;
            end
        end
        observ = select_idx;
    end

    [stageProbs,bundlePI,bundleJ,OmegaInfo,outk] = ...
            naiveSDCore(stageProbs,rpm,OmegaInfo,bundlePI,bundleJ,observ,k,constL,opts,x1);
    x1 = [];
    out.lb(k,1) = outk.lb;
    out.x = outk.x;
    
    if mod(k,50) == 0
        fprintf("[SD][It%d] lb = %f.\n",k,outk.lb);
        for t = 2:length(stageProbs)
            probvec = OmegaInfo.weights / OmegaInfo.cnt;
            [~,I] = sort(OmegaInfo.indices);
            if length(rpm.probs) == length(probvec(I))
                KLdist = sum(rpm.probs.*log(rpm.probs./probvec(I)));
            else
                KLdist = inf;
            end
            fprintf("[SD][It%d][stage-%d] KL divergence = %.4e\n",k,t,KLdist);
        end
    end
end

end