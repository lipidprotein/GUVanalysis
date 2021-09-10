% Virtual GUV creation master script
% 
% createGUVimages(_binary) creates groundtruth images
% cgn_run_images coverts groundtruth into simulated confocal images. 
%
% Il-Hyung Lee (leei@montclair.edu), 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mkdir('Virtual/C1');
mkdir('Virtual/C2');

%% Parameters 

Npts=[504 504 504]; %dimension of the ground truth

confpix=[2 2 8]; %coversion ratio between ground trugh and virtual image

%confpsf=[3.25 3.25 5.5];
%confpsf=[3 3 5];
confpsf=[4 4 7]; %dimension of point spread function in pixels



for index = 1:1:999
    disp(['Working on GUV number :' num2str(index)]);
    
    R = 170 + rand() * 60; %radius of the ground truth vesicle
    L = 1.0 + rand(); %bilayer thickness of the ground truth
    
    if rand() < 0.5
        value = 1000 + 1000*rand();
        value2 = round(value* (0.5 * rand()));
    else
        value = 1000 + 1000*rand();
        value2 = 20 * rand();
    end

    noise = 'SampleR.tif';
     
    
    if index<10
        outName = ['Virtual/C1/' 'phase_' '00' num2str(index) '.tif'];
    elseif index<100
        outName = ['Virtual/C1/' 'phase_' '0' num2str(index) '.tif'];
    else
        outName = ['Virtual/C1/' 'phase_' num2str(index) '.tif'];
    end
    
[guvTruth, imageInfo] = createGUVimages_binary(Npts, R, L, value, value2);

cgn_run_images(guvTruth, confpix, confpsf, noise, outName, imageInfo);

end

for index = 1:1:999
    disp(['Working on GUV number :' num2str(index)]);
    R = 170 + rand() * 60; %radius of the ground truth vesicle
    L = 1.0 + rand(); %bilayer thickness of the ground truth
    
    value = 1000 + 1000*rand();
    value2 = value; % never used
    noise = 'SampleR.tif';
    
    
    if index<10
        outName = ['Virtual/C2/' 'hom_' '00' num2str(index) '.tif'];
    elseif index<100
        outName = ['Virtual/C2/' 'hom_' '0' num2str(index) '.tif'];
    else
        outName = ['Virtual/C2/' 'hom_' num2str(index) '.tif'];
    end
    
[guvTruth, imageInfo] = createGUVimages(Npts, R, L, value);

cgn_run_images(guvTruth, confpix, confpsf, noise, outName, imageInfo);
end
