% rng test

load myrng.mat

close all;
clearvars -except notrand
clc; 

rng(notrand)

for i=1:10
    x(i)=rand
end


