% Run this script to read the length in um per pixel

PathName = './Sample/A';
filelistA = dir([PathName, '/*.tif']);
fnameA = [PathName '/' filelistA(1).name];
distance = pixel_read(fnameA);

