clc;
clear all;
fc=1;
% fs=0.25/4.4717568;
fs=1/64;
dt=[0:1/64:2-1/64];
a=1.2*sin(2*pi*fc*dt)+1*sin(2*pi*3*fc*dt)+2*sin(2*pi*5*fc*dt);%+1*cos(2*pi*3*fc*dt)-1*cos(2*pi*4*fc*dt);
b=0.9*sin(2*pi*fc*(0.6)*dt)+1.6*sin(2*pi*3*fc*(0.6)*dt)+0.5*sin(2*pi*5*fc*dt);%+1*sin(2*pi*3*fc*(1+1)*dt)+0.5*sin(2*pi*4*fc*(1+1)*dt);
plot(a);
% hold on;
% plot(b);
% xlabel("时间","FontName","宋体","FontSize",16);
% ylabel("幅度","FontName","宋体","FontSize",16);
hold on;
% c=interp1([0:length(b)-1],b,[0:length(b)*1.0526-1]/(1+0.05),'spline');%
%     b=interp1([1:length(b)],b,[1:length(b)]/(1+0.05),'spline');
% plot(a);
hold on;
% c=resample(b,19,20);
%     figure();
% plot(c,'b+');
%     d=interp1([1:length(a)],a,(1:855/956:length(a)),'spline');%
%     hold on;
%     figure();
%     plot(d,'r.');