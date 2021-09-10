% Modified from the original program for the GUV analysis purposes 
% 2021 Lee

function [] = cgn_run_images(input_guvTruth, input_confpix, input_confpsf, input_noise, input_outName, input_imageInfo)

% to use ConfocalGN to simulate confocal imaging
% Serge Dmitrieff, Nédélec Lab, EMBL 2015-2016
% http://biophysics.fr

%% Ground truth
truth = input_guvTruth;
% Ground truth can be an image file or a stack, in high resolution

save=1;
display=0;
plot=0;

% pix : number of pixels of the ground truth in a confocal voxel
% ground truth is high res compared to confocal
conf.pix=input_confpix;

% Property of the microscope PSF
% Standard deviation of the PSF in the 3 dimensions provided in units of ground truth pixel size
conf.psf=input_confpsf;


%% Here we present two ways of using ConfocalGN
if 1
    % Reading Noise and Signal parameters from image
    sample = input_noise;
    %% Generating the stack from the image
    % Includes pixel noise
    [stacks,offset,achieved_sig,achieved_noise,im,mean_sig,noise]=confocal_generator(truth,conf,sample);
else
    % Using noise and signal values
    % pixel noise from camera
    noise=[556 1.0371e+04 1.0926e+06]';
    % Estimated mean pix value in signal
    mean_sig=1.0022e+03;
    [stacks,offset,achieved_sig,achieved_noise,im]=stack_generator(truth,conf,noise,mean_sig);
end

%% Save simulated image
if input_imageInfo(1) == 0 && input_imageInfo(2) == 0
    save = 0;
else
    save = 1;
end
if save
    for i =1:1:(floor(input_imageInfo(2)/input_confpix(3)) - ceil(input_imageInfo(1)/input_confpix(3)) + 1)
        temp = ceil(input_imageInfo(1)/input_confpix(3))+i-1;
        % Use if block to remove images at certain z-cooridnates
        %if (temp < 0.15 * size(input_guvTruth,3)/input_confpix(3))||(temp > 0.85 * size(input_guvTruth,3)/input_confpix(3))
        if (temp < 0.20 * size(input_guvTruth,3)/input_confpix(3))||(temp > 0.80 * size(input_guvTruth,3)/input_confpix(3))
            % do not save defined by this range (to avoid images at too low
            % or high parts
        else
            output = strrep(input_outName, '.tif', ['_' num2str(i) '.tif']);
            saveastiff(stacks(:,:,temp),output);
        end
        
    end
end

saveShallow = 0; %save vague images at the top or bottom of GUVs as V_ ...
if saveShallow
    for i=1:1: ceil(size(input_guvTruth,3)/input_confpix(3))
        if i < 0.15 * size(input_guvTruth,3)/input_confpix(3) ||i > 0.85 * size(input_guvTruth,3)/input_confpix(3)
        output = ['V_' strrep(input_outName, '.tif', ['_' num2str(i) '.tif'])];
        saveastiff(stacks(:,:,i),output);
        end
    end
end

%% Display simulated image:
if display
    for i = 1:size(stacks, 3)
        show(stacks(:,:,i));
    end
end

%% Display Image segmentation & plotting to compare simulated data to original
if plot
    plot_simul_results(truth,im,conf,offset)
end
disp('Confocal image generation completed..');
end
