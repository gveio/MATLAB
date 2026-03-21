function [G, H] = getGH_BCH(n, k)
    P = [];
    for i = 1:k
        w = zeros(k, 1);
        w(i) = 1;
        p = bchenc(gf(w'), n, k);
        P = [P; p.x(k+1:n)];
    end
    P = double(P);
    G = [eye(k), P];
    H = [P', eye(n-k)];
end