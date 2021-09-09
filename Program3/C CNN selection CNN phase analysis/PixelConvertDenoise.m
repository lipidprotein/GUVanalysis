%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image converter
% This program converts target images into desired pixels and perform
% thresholding (triangle
% Eli Il-Hyung Lee, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
pixel = 50
    
        
        Path = ['./Virtual/C1/'];
        disp(['Working on foler: ', Path]);
        mkdir('Training/C1');
        tifList = dir([Path, '/*.tif']);
        
        for index = 1:1:length(tifList)
            
            index 
            
            image = imread([Path,  tifList(index).name]);
            
            %data type conversion for compatibility
            image = uint16(image);
                [lehisto x] = imhist(image,11107); %triangle thresholding  
                maskThres = 65535 * triangle_th(lehisto, 11107);
            
     image = image - maskThres;     
     image(image < 0) = 0;

            
            %%Resize into the desired size
            imageOut = imresize(image,[pixel,pixel]);
            
            imwrite(imageOut, ['./Training/C1/' tifList(index).name] );
            
        end
        
        
        Path = ['./Virtual/C2/'];
        disp(['Working on foler: ', Path]);
        mkdir('Training/C2');
        tifList = dir([Path, '/*.tif']);
        
        for index = 1:1:length(tifList)
            
                       
            image = imread([Path,  tifList(index).name]);
            
            %data type conversion for compatibility
            image = uint16(image);
                [lehisto x] = imhist(image,11107); %triangle thresholding  
                maskThres = 65535 * triangle_th(lehisto, 11107);
                
     image = image - maskThres;     
     image(image < 0) = 0;

            
            %%Resize into the desired size
            imageOut = imresize(image,[pixel,pixel]);
            
            imwrite(imageOut, ['./Training/C2/' tifList(index).name] );
            
        end
        
        
        

        
        
        
         