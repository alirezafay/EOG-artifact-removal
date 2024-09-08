clc; close all; clear all; 
sig = load('Ex1_data.mat');
x_org = sig.X_org;
F = abs(fft(x_org));
C = cov(transpose(F));
f1 = 10*10000/100; % fk = (K/N)*fs
f2 = 15*10000/100;
S_x = F(:,f1:f2);
Px = cov(transpose(S_x));
Pxx = 1/2 * (Px + transpose(Px));
[V , D] = eig(Pxx,C);
w1 = V(:,8);
for i=1:8
    V1(:,i)=V(:,9-i);
end
S3 = transpose(w1) * x_org;
S3_den = zeros(8,10000);
S3_den(1,:)= S3;
x3_den = inv(transpose(V1)) * S3_den;
RM = RRMSE(x3_den,sig.X3);
load('Electrodes') ;
offset = max(abs(S3(:))) ;
feq = 100;
ElecName = Electrodes.labels ;
disp_eeg(S3,offset,feq,ElecName);
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
