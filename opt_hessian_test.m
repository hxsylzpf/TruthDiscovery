clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. obj function and grad are provided in bigtoleft(x)
% 2. constraints and grad are provided in twocone(x)
% 3. hessian of obj func and constraints are provided in hessinterior(x)

   % Set options as follows:

   % grad constraint (nonlcon) is used
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
            'Display','off','GradObj','on','GradConstr','on',...
            'Hessian','user-supplied','HessFcn',@hessinterior);
        
    % Run fmincon with starting point [¨C1,¨C1,¨C1], using the options structure:

    tic
    [x fval mflag output] = fmincon(@bigtoleft,[-1,-1,-1],[],[],[],[],[],[],@twocone,options)
    time = toc
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. obj function, grad and hessian are provided in one function.
options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','GradObj','on','Hessian','user-supplied','HessFcn',@rosen_hessian);

% Run fminunc starting at [-1;2]:
lb = [0; 0];
ub = [1;1];

[x, fval] = fmincon(@rosen_grad,[-1,2],[],[],[],[],lb,ub,[],options)
