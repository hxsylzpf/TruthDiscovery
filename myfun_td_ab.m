function f = myfun_td_ab(x, inst, g, z_post)

n = length(inst);

f = 0;
for j = 1:n
        f = f + z_post(j)*(inst(j)*mylog(g(j)*x) + (1-inst(j))*mylog(1- g(j)*x));
end

f = -f;
