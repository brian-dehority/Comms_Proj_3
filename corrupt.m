function corrupted = corrupt(M,p)
    random = rand(size(M));
    random(random < p) = 0;
    random(random >= p) = 1;
    corrupted = mod(M + random,2);
end

   