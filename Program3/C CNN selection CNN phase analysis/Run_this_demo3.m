%%
% Program 3 
% CNN based vesicle selection + CNN based classification of final phase domain separation.
% Phase domain clssification training set is generated completely by virtual
% simulation + training set augmentation. 
%
% This code was written for MATLAB R2020b.
%
% Running this script will show an example analysis using program 3
% In this example, single color stacks with assorted phase separation
% images will be analyzed (fluorophore is TopFluor-Cholesterol)
% Pre-trained neural network information is included in the file 'nrFilter4C.mat'
% which is used for vesicle selection filtering. 
% This demo includes creation of virtual training set and training a fresh
% CNN to use in the example image classification analysis. 
% Due to the fact that training set generation and CNN training relies on
% random generation, the result will not be entierly identical each run, but
% it is most likely show similar capability of analysis at the end. 
% The training stage can be skipped to use the pre-trained network. 
%
% Il-Hyung Lee (leei@montclair.edu), 2021
%%

tic;
clear all;

owd = pwd;
addpath([owd '/VIrtual GUV creator/'])
addpath([owd '/VIrtual GUV creator/cgn/'])
addpath([owd '/VIrtual GUV creator/matned/'])
addpath([owd '/VIrtual GUV creator/segmentation/'])
addpath([owd '/VIrtual GUV creator/GUV_ground_truth_generator/'])

%% Creating virtual GUV training set. (can be skipped and use provided CNN)
%% Running this will take some time and will create >10,000 files 

run_GUV;
PixelConvertDenoise;
clear all;
Train;

%%

%% This part actually performs vesicle analysis 
C_NeuralPhase;
%%

clear all;
toc;