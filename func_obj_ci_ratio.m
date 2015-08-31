function f = func_obj_ci_ratio(x, n_pt, n_s, X)

% x is constructed as [z a b]
% z - x(1) to x(n_s)
% a - x(n_s + 1) to x(n_s + n_pt)
% b - x(n_s + n_pt + 1) to x(n_s + 2*n_pt)

f = 0;
for j = 1:n_s
    for i = 1:n_pt
        if ~isnan(X(i,j))
            f = f + x(n_s+i)*mylog(x(n_s + n_pt + i)) + (x(n_s+i) - 1)*(mylog(X(i,j)) - mylog(x(j))) - x(n_s + n_pt + i)*X(i,j)/x(j);
        end
    end
end

f = -f;

% if nargout > 1
%     grad_z = zeros(n_s, 1);
%     grad_lambda = zeros(n_pt, 1);
%     
%     for j = 1:n_s
%         grad_z(j) = - prior_nu(j)*(x(j) - prior_mu(j));
%         for i = 1:n_pt
%             if X(i,j)~=0
%                 grad_z(j) = grad_z(j) - x(n_s+i)*(x(j) - X(i,j));
%             end
%         end
%     end
%     
%     for i = 1:n_pt
%         grad_lambda(i) = (prior_a(i) - 1)/x(n_s+i) - prior_b(i);
%         for j = 1:n_s
%             if X(i,j)~=0
%                 grad_lambda(i) = grad_lambda(i) + 0.5/x(n_s+i) - 0.5*(X(i,j) - x(j))^2;
%             end
%         end
%     end
%     
%     grad = [grad_z; grad_lambda];
%     grad = -grad;
% end

