clc; close all; clear all; 
sig = load('pure.mat');
sig2 = load('contaminated.mat');
pure1 = sig.pure;
contaminated = sig2.contaminated;
time_EOG = [1.5 2.925; 16 16.5; 21.9 22.5];
tt = time_EOG *200; % for defining sample numbers 
x_hat = zeros(19,508);
C = cov(transpose(pure1));
x_hat(:,1:286)=contaminated(:,300:585); % defining EOG activity Signal
x_hat(:,287:387)=contaminated(:,3200:3300);
x_hat(:,388:508)=contaminated(:,4380:4500);
C_hat = cov(transpose(x_hat));
[V , D] = eig(C_hat,C);
w1 = V(:,1:17); % defining Vectors that are not related to EOG activities
S = transpose(w1) * contaminated;
S_denoised = zeros(19,5601);
S_denoised(1:17,:) = S;
x_den = inv(transpose(V)) * S_denoised;
load('Electrodes') ;
offset = max(abs(x_den(:))) ;
feq = 200;
ElecName =  [] ;
titre = 'denoised signal'; 
disp_eeg(x_den,offset,feq,ElecName,titre);
RM = RRMSE(pure1,x_den);
function [R] = RRMSE(x_org,x_den)
    sum1 = 0;
    sum2 = 0;
    for i=1:length(x_org(:,1))
        a = (x_org(i,:)-x_den(i,:)).^2;
        b = x_org(i,:).^2;
        sum1 = sum1 + sum(a);
        sum2 = sum2 + sum(b);
    end
    R = sqrt(sum1/sum2);
end


