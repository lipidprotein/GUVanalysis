function [outputArg] = pixel_read(filePath_in)
%PIXEL_READ This function reads a tiff metadata to return the value of
%distance per pixel
t = Tiff(filePath_in,'r');
tagValue = t.getTag('XResolution');


outputArg = 1/tagValue;

disp([num2str(outputArg) 'um per pixel']);

end

