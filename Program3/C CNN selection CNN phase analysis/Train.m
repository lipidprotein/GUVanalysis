%%
% train CNN using training images as defined by the structure below
%%

imds = imageDatastore('./Training/', ...  
    'IncludeSubfolders',true,'LabelSource','foldernames');

temp1 = dir('./Training/C1/');
temp2 = dir('./Training/C2/');
temp = min(length(temp1),length(temp2)) -2 % number of files  
numTrainFiles = temp;
[imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,'randomize');

%% augmentation part
inputSize = [50, 50, 1];
pixelRange = [-10 10];
scaleRange = [0.9 1.1];
imageAugmenter = imageDataAugmenter( ...
    'RandXReflection',true, ...
    'RandXTranslation',pixelRange, ...
    'RandYTranslation',pixelRange, ...
    'RandXScale',scaleRange, ...
    'RandYScale',scaleRange);
augimdsTrain = augmentedImageDatastore(inputSize(1:2),imdsTrain, ...
    'DataAugmentation',imageAugmenter);
%%


layers = [
    imageInputLayer([50 50 1])
    
    convolution2dLayer(3,8,'Padding','same')
    %convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    %convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding','same')
    %convolution2dLayer(3,64,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',50, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',30, ...
    'Verbose',false, ...
    'Plots','training-progress');


%net2 = trainNetwork(imdsTrain,layers,options);
net2 = trainNetwork(augimdsTrain,layers,options);


YPred = classify(net2,imdsValidation);
YValidation = imdsValidation.Labels;

accuracy = sum(YPred == YValidation)/numel(YValidation)

%save neural network (by saving variables) 
save nrFilter2Cstate_net2;
