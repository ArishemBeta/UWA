clear all;
close all;
clc;

L=1000;
parpool('local',4);
tic
X=rand(64,64,L);
Y=zeros(64,64,L);
parfor i= 1:L
    i
    for j=1:200
        Y(:,:,i) = inv(X(:,:,i));
    end
end
toc
delete(gcp('nocreate'));