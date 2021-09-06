%%
% Program 1 
% All computation based binding intensity analysis. 
%
% This code was written for GNUOctave 6.3.0.
%
% Running this script will show an example analysis using program 1
% In this example, two color stacks will be analyzed 
% Channel A is reference and Channel B is for protein bound to the membrane
% One stack set has high ChB intensity (high binding) 
% and another stack set has low ChB intensity (low binding) for comparison
% detailed results, if needed will be saved in a folder 'Sample/Analysis'
%
% Il-Hyung Lee (leei@montclair.edu), 2021
%%

tic;
clear all;
pkg load image;

A_Intensity

clear all;
toc;