function ff=switchingtime(a,b)
% generate switching time for rate of the form max(0, a + b s) + c
% under the assumptions that b > 0, c > 0,
u=rand;
if (b > 0)
    if (a < 0)
        ff= -a/b + switchingtime(0.0, b);
    else % a >= 0
        ff=-a/b + sqrt(a^2/b^2 - 2 * log(u)/b);
    end
elseif (b == 0) % degenerate case
    if (a < 0)
        ff= Inf;
    else % a >= 0
        ff= -log(u)/a;
    end
else % b <= 0
    if (a <= 0)
        ff= Inf;
    else % a > 0
        y = -log(u); t1=-a/b;
        if (y >= a * t1 + b *t1^2/2)
            ff= Inf;
        else
            ff= -a/b - sqrt(a^2/b^2 + 2 * y /b);
        end
    end
end
end
