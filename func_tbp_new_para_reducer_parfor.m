function record = func_tbp_new_para_reducer_parfor(mapper_files, n_node, n_pt, n_l, prior_a, prior_b)

n_s = 0;
s0 = 0;
s1 = 0;
s2 = 0;
s3 = 0;
z_est = [];

% reduce
for m = 1:n_node
    cur_record = load(mapper_files{m});
    n_s = n_s + cur_record.n_s;
    s0 = s0 + cur_record.s0;
    s1 = s1 + cur_record.s1;
    s2 = s2 + cur_record.s2;
    s3 = s3 + cur_record.s3;
    z_est = [z_est cur_record.z_est];
end

% derive
pi_est = (s0 + 1)/(n_s + n_l);

h_est = zeros(n_pt, n_l);
lambda_est = zeros(n_pt, n_l);

for i = 1:n_pt
    h_est(i,:) = s2(i,:)./s1(i,:);
    lambda_est(i,:) = (0.5*s1(i,:) + prior_a(i,:) - 1)./(0.5*s1(i,:).*h_est(i,:).^2 - s2(i,:).*h_est(i,:) + 0.5*s3(i,:) + prior_b(i,:));
end

record.pi_est = pi_est;
record.h_est = h_est;
record.lambda_est = lambda_est;
record.z_est = z_est;

