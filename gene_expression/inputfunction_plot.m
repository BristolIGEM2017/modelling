% Script to plot steady-state promoter activity
% Jeremie Joannes
% 31/08/2017
close all; clear all; clc;

LacI = 1e-6; % (Alon)

IPTG = linspace(0,1,1000);

B = 1e-3; % mRNA/s max transcription rate

Kx = 1e-6; % M (Alon)
Kx_on = 1e9; % 1/M/s (Alon)
Kx_off = Kx*Kx_on; % 1/s

Kd = 1e-10; % M (http://kirschner.med.harvard.edu/files/bionumbers/Variabl
% es%20and%20constants%20used%20for%20modeling%20the%20lac%20operon.pdf)
Kd_on = 1e9; % 1/M/s (Alon)
Kd_off = Kd*Kd_on; % 1/s


f = B./(1+LacI./Kd./(1+IPTG./Kx));

figure(1)
plot(IPTG,f)
title('Steady-state input function')
xlabel('[IPTG] in M')
ylabel('Promoter activity in mRNA/s')

