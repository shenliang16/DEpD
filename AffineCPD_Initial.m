opt.whatHapen=0;
if ~exist('P_thr','var') || isempty(P_thr) || (P_thr==0);  P_thr=0.75;  end
if ~isfield(opt,'max_it') || isempty(opt.max_it), opt.max_it = 500; end;
if ~isfield(opt,'tol') || isempty(opt.tol), opt.tol = 5e-4; end;
if ~isfield(opt,'viz') || isempty(opt.viz), opt.viz = 0; end;
if ~isfield(opt,'corresp') || isempty(opt.corresp), opt.corresp = 0; end;
if ~isfield(opt,'outliers') || isempty(opt.outliers), opt.outliers = 0.5; end;
if ~isfield(opt,'fgt') || isempty(opt.fgt), opt.fgt = 0; end;
if ~isfield(opt,'sigma2') || isempty(opt.sigma2), opt.sigma2 = 0; end;

% strictly rigid params
if ~isfield(opt,'rot') || isempty(opt.rot), opt.rot = 1; end;
if ~isfield(opt,'scale') || isempty(opt.scale), opt.scale = 1; end;
% strictly non-rigid params
if ~isfield(opt,'beta') || isempty(opt.beta), opt.beta = 3; end;
if ~isfield(opt,'lambda') || isempty(opt.lambda), opt.lambda = 3; end;
% lowrank app param
if ~isfield(opt,'numeig') || isempty(opt.numeig), opt.numeig = round(sqrt(N)); end;
if ~isfield(opt,'eigfgt') || isempty(opt.eigfgt), opt.eigfgt = 1; end;

