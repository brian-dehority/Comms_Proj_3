function corrected = correct(map, H, rec);

corrected = zeros(length(rec));
for i = 1:length(rec)
    k = rec(i,:); %single received word
    mult = f2mult(k, H.'); %multiply by parity check matrix
    error = map(bi2de(mult)); %convert to decimal, to find it in the map
    corrected(i,:) = f2add(k,de2bi(error, length(k))); %correct word
end 
end
