function corrected = correct(map, H, rec);
corrected = zeros(size(rec)); 

for i = 1:size(rec,1)
    k = rec(i,:); %single received word
    
    mult = f2mult(k, H.'); %multiply by parity check matrix
    multbi = bi2de(mult);
    if isKey(map,multbi) == 1 %Checks to see if the syndrome is in the map
        error = map(multbi); %convert to decimal, to find it in the map
        corrected(i,:) = f2add(k,de2bi(error, length(k))); %correct word
    else %if syndrome not in map
        corrected(i,:) = NaN * ones(size(k)); %fill with NaN
    end
end 
end
