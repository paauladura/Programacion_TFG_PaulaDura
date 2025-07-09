function ordrng = itsg_sh_ordering(lmax, start_deg)
    
ordrng  = zeros((lmax+1)^2-4, 2);

n = 0;
for l = start_deg:lmax
    for m = 0:l
        n = n + 1;
        if m == 0
            ordrng(n, :) = [l, m];
        else
            ordrng(n, :) = [l, m];
            
            n = n + 1;
            ordrng(n, :) = [l, -m];
        end
    end
end
