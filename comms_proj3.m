%% generates the syndrome arrays 
clc; clear; close all;

%The Generator and Code book matricies

G1 = [1, 0, 0, 0, 0, 1, 1, 1;
      0, 1, 0, 0, 1, 1, 1, 0;
      0, 0, 1, 0, 1, 0, 1, 1;
      0, 0, 0, 1, 1, 1, 1, 1;]
  
G2 = [1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0;
      0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0;
      0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1;
      0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1]
  
C1 = [0,0,0,0,0,0,0,0;
     0,0,0,1,1,1,1,1;
     0,0,1,0,1,0,1,1;
     0,0,1,1,0,1,0,0;
     0,1,0,0,1,1,1,0;
     0,1,0,1,0,0,0,1;
     0,1,1,0,0,1,0,1;
     0,1,1,1,1,0,1,0;
     1,0,0,0,0,1,1,1;
     1,0,0,1,1,0,0,0;
     1,0,1,0,1,1,0,0;
     1,0,1,1,0,0,1,1;
     1,1,0,0,1,0,0,1;
     1,1,0,1,0,1,1,0;
     1,1,1,0,0,0,1,0;
     1,1,1,1,1,1,0,1]
 
C2 = [0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,1,1,0,1,0,1,1,1,1;
      0,0,1,0,1,1,0,1,1,0,1,1;
      0,0,1,1,0,1,1,1,0,1,0,0;
      0,1,0,0,1,0,0,1,1,1,1,0;
      0,1,0,1,0,0,1,1,0,0,0,1;
      0,1,1,0,0,1,0,0,0,1,0,1;
      0,1,1,1,1,1,1,0,1,0,1,0;
      1,0,0,0,1,1,1,1,0,1,1,0;
      1,0,0,1,0,1,0,1,1,0,0,1;
      1,0,1,0,0,0,1,0,1,1,0,1;
      1,0,1,1,1,0,0,0,0,0,1,0;
      1,1,0,0,0,1,1,0,1,0,0,0;
      1,1,0,1,1,1,0,0,0,1,1,1;
      1,1,1,0,1,0,1,1,0,0,1,1;
      1,1,1,1,0,0,0,1,1,1,0,0]

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
hamming1 = size(G1, 2)*min(pdist(G1, 'hamming'))
%This is one way of doing hamming, but idk if its right ?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% code 1
%defining constants
k1 = 4;
n1 = 8;
stand_1 = zeros(n1+1, 2^k1, n1);

stand_1(1, :, :) = C1;

I = eye(n1);

%Standard errors
for i = 1:n1
    stand_1(i+1, :, :) = f2add(C1, I(i, :));
end

%Parity matricies
[kg1, ng1] = size(G1);
parityBits1 = G1(:, (kg1 + 1:end));
H1 = [transpose(parityBits1) eye(ng1-kg1)];
    
   
% creating the syndrome and error matricies
E1 = reshape(stand_1(:, 1, :), [], n1);
S1 = f2mult(E1, H1.');
E1dec = bi2de(E1);
S1dec = bi2de(S1);
%The map seems like the easiest way to relate the syndromes to errors
StandE_1 = containers.Map(bi2de(S1), bi2de(E1)); %we can store the vlaues as decimal in a map

%% code 2
%Same process as above
k2 = 4;
n2 = 12;
stand_2 = zeros(1+12+12*11/2, 2^k2, n2);

% 12
I = eye(n2);

stand_2(1, :, :) = C2;

counter = 2;
for i = 1:n2
    stand_2(counter, :, :) = f2add(C2, I(i, :));
    counter = counter + 1;
end

for i = 2:n2
    for j = 1:i-1
        stand_2(counter, :, :) = f2add(C2, I(i, :) + I(j, :));
        counter = counter + 1;
    end
end        

%Parity matricies
[kg2, ng2] = size(G2);
parityBits2 = G2(:, (kg2 + 1:end));
H2 = [transpose(parityBits2) eye(ng2-kg2)];
 




E2 = reshape(stand_2(:, 1, :), [], n2);
S2 = f2mult(E2, H2.');
E2dec = bi2de(E2);
S2dec = bi2de(S2);
StandE_2 = containers.Map(S2dec, E2dec); %Creates a map between the syndromes and errors for G2
%%


