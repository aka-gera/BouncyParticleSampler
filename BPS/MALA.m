function Xi=MALA( Xi_init, h,  n_iter,fun)

% % grad_U=@(x)fun{2,3}(x);
% % U_fun=@(x)fun{1,3}(x);
Xi = zeros(2,  1);
Xi(:,1) = Xi_init;
j1=3;
q=@(x, y)  exp(-1/(4*h) * norm((y -( x - (h)*fun{2,j1}(x))))^2);
acceptance_count = 0;
tic;
for i=1:n_iter
    Xi_current = Xi(:,i);
    Xi_proposed = Xi_current -(h) *fun{2,j1}(Xi_current) + sqrt(2*h) * randn(2, 1);
    accept_prob = min(1, exp(fun{1,3}(Xi_current)-fun{1,3}(Xi_proposed)) * q(Xi_proposed, Xi_current) /q(Xi_current,Xi_proposed));
    if rand <= accept_prob
        Xi(:,i+1) = Xi_proposed;
        acceptance_count =acceptance_count+ 1;
    else
        Xi(:,i+1) = Xi_current;
    end
end
ttime=toc;
disp(['MALA: Fraction of accepted proposals: ', num2str(acceptance_count / n_iter)])
disp(['MALA: Time of simulation: ', num2str(ttime)])
% dlmwrite('xMala.txt',Xi','delimiter','\t','precision',5)
end