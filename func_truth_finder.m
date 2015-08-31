function record = func_truth_finder(n_state, n_pt, n_event, provided_label_mat, event_label, loc_cover_mat)

n_iter = 100;

gamma = 0.3;
eps = 1e-4;

% initialization
t0 = 0.8*ones(n_pt, 1);

% number of locations covered by a pt
% F = sum(loc_cover_mat, 2);
% F = sum(provided_label_mat, 2);
F = n_event*ones(n_pt, 1);

A = zeros(n_pt, 2*n_event);

B = zeros(2*n_event, n_pt);

for i = 1:n_pt
    for j = 1:n_event
       % provide a fact
       if loc_cover_mat(i,j) == 1
           if provided_label_mat(i,j) == 0
               A(i, 2*j-1) = 1/F(i);
               B(2*j-1, i) = 1;
           elseif provided_label_mat(i,j) == 1
               A(i, 2*j) = 1/F(i);
               B(2*j, i) = 1;
           end
       end
    end
end

t = t0;

for i = 1:n_iter
    t_c = t;
    tau = - log(1 - t);
    sigma = B*tau;
    % confidence score of facts
    s = 1./(1 + exp(- gamma*sigma));
    % t must be smaller than 1
    t = A*s;
    
    if t.'*t_c/(norm(t)*norm(t_c)) >= 1 - eps
        break;
    end
end

% initiated as 0
est_label = zeros(1, n_event);
% generate the final label for event
for i = 1:n_event
    if s(2*i - 1) < s(2*i)
        est_label(i) = 1;
    end
end

% confusion matrix
[~, confmtx] = confmat(n_state, n_event, event_label, [0 1], est_label)
[~, ~, ~, pre1, rec1, F, WAF] = binary_f_measure(confmtx)
   
    record.confmtx = confmtx;
    record.F = F;
    record.WAF = WAF;
    record.est_label = est_label;
    record.pre = pre1;
    record.rec = rec1;
    