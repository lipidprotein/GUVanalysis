% GUV groundtruth generator. 
% 07/16/2021 by Eli Lee
% This program uses user input parameters to create artificial GUVs. 
% This version returns information to selectively save single images that 
% contains enough circular intensity. 

function [guvTruth, imageInfo ] = createGUVimages_binary(input_Npts, input_R, input_L, input_value, input_value2)
%% Variables defined

% Name of the file to save
% fileName = 'output.tiff';

% Size of the ground truth in discretized units
Npts=input_Npts;

% Radius of the GUV in discretized units (shouldn't be greater than the
% image size)
R = input_R;

% Thickness of the lipid membrane in discretized units. Lx2 will be the
% bilayer thickness
L = input_L;

% Intensity value to be asigned for the primary intensity
value = input_value;



%% GUV creation

% Initialization
guvTruth = zeros(Npts);

% Options for complex GUVs
% Center position shifter x-y on/off (to prevent all images exactly
% centered, choosing this will move the center of the vesicle in x-y.
isCenterShift = 1;
% amount of shift, will be determined randomly from +/- zero to this value.
% Avoid Radius + shift being greater than 1/2 of the dimension (in pixels
% of ground truth)
centerShift = 20;
% binary phase separation on/off
isBinary = 1;
% Intensity value to be asigned for the second domain intensity
value2 = input_value2;
% Radius of the secondary domain (larger value means greater area)
binaryRadius = R * 1.0;
% Amount of random variation in the binary Radius (1.0 is 100%, larger
% value means more variations in relative domain sizes) 
binaryVariance = 0.5;

% true speckel noise addition on/off
isNoise = 1;
% how much to add 0.01 means 1% of the pixesl will be added as spekcle noise
noiseRatio = 0.00001; 
% Noise intensity, noise intensity will be randomized from zero to this
% value
noiseIntensity = value; 

% small vesicle addition on/off
isVesicle = 1;
% how many vesicles to add per image stack
vesicleNumber = 100; 
% Vesicle intensity,intensity will be randomized from zero to this value
vesicleIntensity = value; 
% Vesicle radius. will be randomized from zero to this value
vesicleRadius = 30;

% ratio of pixels in one domain should be at least this ratio of the other.
% threshold applied for separate image saving information. 
% cutThres = round(Npts(1)*Npts(2)*0.005 * 0.1);
%cutThres = round(Npts(1)*Npts(2)*0.01*0.4);
cutThres = round(Npts(1)*Npts(2)*0.01*0.2);


%% center coordinate of the GUV created
if isCenterShift
    center = [Npts(1)/2 + 2.0 * (rand()-0.5)*centerShift , Npts(2)/2 + 2.0 * (rand()-0.5)*centerShift, Npts(3)/2];
else
    center = [Npts(1)/2, Npts(2)/2, Npts(3)/2];
end
%%

%% Scan through all voxels and give a value if it is within the definition
% of the GUV membrne (brute force calculation for algorithmic simplicity)
for z = 1:1:Npts(3)
    for y = 1:1:Npts(2)
        for x = 1:1:Npts(1)
            dist = norm( [x,y,z] - center);
            if dist > (R - L)  &&  dist < (R + L)
                guvTruth(x,y,z) = value;
            end
                       
        end
    end
end

disp('Vesicle succefully created..');
%%

%% Create the second phase separation domain based on the setting. Another
% globe on a random point on the surface of the original GUV will be
% created. Radius of this domain will be determined randomly based on the
% given parameters. 

if isBinary
    tempX = 2*(rand()-0.5);
    tempY = 2*(rand()-0.5);
    tempZ = 2*(rand()-0.5);
    tempNorm = norm([tempX, tempY, tempZ]);
    
    binaryCenter(1) = center(1) + tempX * R / tempNorm;
    binaryCenter(2) = center(2) + tempY * R / tempNorm;
    binaryCenter(3) = center(3) + tempZ * R / tempNorm;
    
    tempRadius = (1.0 + 2.0 * (rand()-0.5) * binaryVariance) * binaryRadius;
    
    for z = 1:1:Npts(3)
        for y = 1:1:Npts(2)
            for x = 1:1:Npts(1)
                %only consider when the voxel has the primary intensity value
                if guvTruth(x,y,z) == value
                    dist = norm( [x,y,z] - binaryCenter);
                                if dist < tempRadius
                                    % second globe is much greater to
                                    % reduce chopping effect? 
                                    %circle in polar coordinate
                                    guvTruth(x,y,z) = value2;
                                end
                end
                
            end
        end
    end
end

disp('Vesicle separation completed..');

%%
%% imageInfo identifier - calculating the range of z that contains enough images
% to be used as single training images. 
imageInfo = [0, 0];

for z=1:1:Npts(3)
    numOfPhase1 = sum(sum(guvTruth(:,:,z)==value));
    numOfPhase2 = sum(sum(guvTruth(:,:,z)==value2));
    if (numOfPhase1 > cutThres) && (numOfPhase2 > cutThres)
        if imageInfo(1) == 0
            imageInfo(1) = z;
            imageInfo(2) = z;
        else
            imageInfo(2) = z;
        end
    end
end
%Augmentation by adding variations
%Center shifter 



%% Vesicle addition unit
if isVesicle 
    for i = 1:1:vesicleNumber
        
    tempX = ceil(rand()*Npts(1));
    tempY = ceil(rand()*Npts(2));
    tempZ = ceil(rand()*Npts(3));
    tempInt = ceil(rand()*vesicleIntensity);
    tempRad = ceil(rand()*vesicleRadius); 
    
    for z = round(tempZ - tempRad - L):1:round(tempZ + tempRad + L)
        for y = round(tempY - tempRad - L):1:round(tempY + tempRad + L)
            for x = round(tempX - tempRad - L):1:round(tempX + tempRad + L)
                %only consider when the voxel has the primary intensity value
                if (x<1)||(y<1)||(z<1)||(x>Npts(1))||(y>Npts(2))||(z>Npts(3))
                    continue;
                else
  
                    dist = norm( [x,y,z] - [tempX, tempY, tempZ]);
                                if dist > (tempRad - L)  &&  dist < (tempRad + L)
                                    guvTruth(x,y,z) = guvTruth(x,y,z) + tempInt;
                                end

                end
                
            end
        end
    
    end
    end
end
disp('Vesicle addition ccompleted..');
%%


%% True noise adder unit
% This part add true speckel fluorescence noise to the image based on the
% parameters given. 
if isNoise
for i = 1:1:round(noiseRatio * Npts(1)* Npts(2)*Npts(3))
   
    x = ceil(rand()*Npts(1));
    y = ceil(rand()*Npts(2));
    z = ceil(rand()*Npts(3));

    guvTruth(x,y,z) = guvTruth(x,y,z) + ceil(rand()*noiseIntensity);
    
end
end
disp('Noise speckel addition completed..');
%% 




%%
guvTruth = single(guvTruth);
%saveastiff(guvTruth, fileName);
end

