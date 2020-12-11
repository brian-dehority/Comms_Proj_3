function corrupted = corrupt(M,p)
    N = length(M);
    random = rand(N,1);
    random(random < p) = 0;
    random(random >= p) = 1;
    corrupted = mod(M + random,2);
end

   