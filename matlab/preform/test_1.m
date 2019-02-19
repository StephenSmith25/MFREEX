clear all
close all

sigma = [2 1 0 ; 1 3 0 ; 0 0 4]


S = sigma - (1/3)*trace(sigma)*eye(3);


R  = 0;
for i = 1:3
    for j = 1:3
        R = sigma(i,j)*sigma(i,j);
    end
end
R