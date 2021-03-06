function [ stack,offset] = generate_stacks( img,conf,sig,noise,options)
% generate_stacks : make mock confocal data from a ground truth
%   Distributed under the terms of the GNU Public licence GPL3
%
% INPUT
% * img is the orginal ground truth, with isotropic 1x1x1 pixels
%   this image shall be convolved with the psf
%   should already have background noise included if any
%   is either a matrix or a tiff file
%
% * conf.psf is the point spread function of the microscope 
%       it is the 2-way psf (i.e. illumination+observation)
%       also called confocal psf  ; here : 1x3 vector
%   here we assume the psf to be gaussian which provides fast convolution
%   the psf is in units of the size of original image pixel 
%       ex : psf = [wx,wy,wz] where wi is the gaussian width on axis i
%
% * conf.pix is the voxel size of the microscope, i.e. how many pixels of the
%   original image fit in one microscope voxel         (1x3 vector)
%
% * noise is (at least) the two first moments of the noise distribution
%   If the distribution is skewed enough (large enough third moment)
%   it is modeled by gamma function
%   Otherwise, the noise is gaussian
%
% * sig is the moments of the signal distribution
%       sig(1) is the aim for the signal mean pixel value
%
% OUTPUT
% * stack is the simulated confocal stack
%
% * offset is the translation in x,y,z to go back to the original coord
%   as the shape gets
%
% RATIONALE
% We convolve the ground truth with the confocal psf. The psf should be the
% one measured experimentally from control sample (e.g. bead). So it is
% really how a point *object* would be seen.
%
% In many cases, the confocal psf is close to gaussian. Please make sure it
% is the case for your setup, or that the exact shape of the PSF is not
% relevant to your issue.
%
% This program is intended solely for analysis testing purpose.
% This IS not a rigorous implementation of confocal optics and electronics.
%
% Serge Dmitrieff, Nédélec Lab, EMBL 2016
% www.biophysics.fr

if nargin<2
    error('You must provide an image and confocal imaging properties');
end
psf=conf.psf;
pix=conf.pix;

% If img is the name of a tiff file, importing it as an array
if ischar(img)
    imgs=tiffread(img);
    img=get_stack(imgs);
end

if nargin>2
   if sig(1)<0
        error('Invalid expected mean signal');
   elseif sig(1)<noise(1);
       warning('Signal expected lower than signal');
   end
else 
    error('No signal level specified');
end


if nargin<4
    noise=[0 0 0];
    warning('No noise');
else
    if length(noise)<2
        error('Invalid noise properties (not enough moments)');
    else
        if length(noise)==2
            noise(3)=0;
        end
    end
end


if nargin<5
    options=cgn_options_load();
end
if isfield(options,'segmentation')
    options=options.segmentation;
end


%% 3D Convoluting of the image
% by default, convolved with a gaussian kernel

img=convolve_with_psf(img,psf);

%% Taking the pixels 

[stack,offset,nn]=stack_from_img(img,pix);

%% Generating pixel noise
px_noise=pixel_noise(noise,nn);

%% Estimating signal & noise levels and calibrating noise
[ Ssig,Snoise] = get_img_params(stack,options);
f=Ssig(1);
b=Snoise(1);
% Desired signal to noise 
s=sig(1)/noise(1);
% Linear algebra on stack to reach desired signal and noise level
B=(f-s*b)/(s-1);
A=1/(b+B);
stack=A*(stack+B);
%% We multiply the stack by the noise
stack(:)=stack(:).*px_noise(:);
end

