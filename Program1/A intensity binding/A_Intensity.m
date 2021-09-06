%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUV binding intensity Analysis program 
% This program uses user defined variables (size of vesicles, sensitivity)
% to automatically pick up GUV, and analyze fluoresence intensity around it.
% Multi channel analysis to normalize intensity. 
% In the garget folder, folder named "A" should have channel A images and 
% another folder named "B" should have matching channel B images
% Images should be exact matches and should be sorted alphabetically 
%
% Il-Hyung Lee (leei@montclair.edu), 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instruction: This Code was written for GNUOctave. Running it in Matlab 
% requires modification of the code.
% Target files should be separated into Subfolders of A and B.
% Users can choose the target folder contianing the subfolders. triangle_th.m 
% should be accessible for this program. Images should
% be in Tiff format and A and B should be matching images of different
% color channels. Output files will be automatically generated. 
% Output images show initial circle detection for each z-section and then
% final GUVs analyzed in the final two images for eacn z-stack. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output txt format 
%% Each raw is for one GUV detected
%% from left to right:
%% 1. Average coordinate x of the GUV
%% 2. Average coordinate y of the GUV
%% 3. Maximum radius of the GUV (Greatest raadius of multiple Z-images)
%% 4. Total number of Z-images of the GUV
%% 5. Total sum of intensity near membrane ChA (average of multiple Z-images)
%% 6. Total number of pixels counted near membrane ChA (average of multiple Z-images)
%% 7. Average intensity per pixel ChA (average of multiple Z-images)
% 8. Average background intensity ChA
% 9. Average intnsity per pixel - Average background intensity ChA (7 - 8)
%% 10. Total sum of intensity near membrane ChB (average of multiple Z-images)
%% 11. Total number of pixels counted near membrane ChB (average of multiple Z-images)
%% 12. Average intensity per pixel ChB (average of multiple Z-images)
% 13. Average background intensity ChB
% 14. Average intnsity per pixel - Average background intensity ChB (10 - 11)
%% 15. Normalized ChB intentity = Adjusted Average intentsity B / Adjusted Average intentisy A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important parameters to set by users
radius_min = 20; %smallest radius of GUV to detect
radius_max = 100; %largest radius of GUV to detect
sensitivity = 0.90; %CHT fit sensitivity (0-1, 1 means free pass)
edge = 0.1; %CHT edge definition sensitivity (0-1, 1 means trong edge needed)
thres = 4; %number of pixels to count is within the fluorescence
thres2 = thres^2; 
trackDist = 5; %pixel distance to consider the same vesicles
initSkip = 1; %%number of images to skip (in case initial few are not good)

autoBackground = 1; %automatic estimation of background intensity by most probable 
                    %vlaue 1 to use mode estimation, value 2 for median estimation
backgroundIntA = 0;  %manual input of the backgroudn intensity in case not using auto estimation
backgroundIntB = 0;

saveCircles = 1; %1 - on to save all circle figures as jpeg 
% Important parameters to set for GUV stack analysis 
autoBinaryMask = 1; % 1 to use automatic binary mask thesholding 2 to use manual mask value 0 to use raw images
maskThres = 0.5; % if not using automask, any intensity greater than this avlue will be 1 and smaller will be 0

isMLVfilt = 1; % 1 for filtering MLVs 0 for no filer
MLVthres = 2.0; % 1.5 means outer membrane adjusted intensity should be at least 150% inner membrane adjsuted intensity to pass
minNumStack = 5; % Minimum number of z-stacks needed to count as GUV (0 for not using it)

isCloseStackfilt = 2; % 1 for filtering GUV Stacks that are too close
                      % 1 to remove both GUVs stuck 2 to reject only one with less fitting fitness

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PathName=uigetdir('Choose file folers containing directory (two colors)');
PathName = './Sample/';
filelistA = dir([PathName, '/A/*.tif']); % reference channel
filelistB = dir([PathName, '/B/*.tif']); % binding or changing channel
if length(filelistA) != length(filelistB)
  disp('number of stacks mismatch in two channels');
  exit;
endif
PathOut = [PathName '\Analysis'];
if exist([PathName '\Analysis']) != 7
  mkdir(PathOut);
endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=1:1:length(filelistA)
 
  fnameA = [PathName '\A\' filelistA(index).name];
  fnameB = [PathName '\B\' filelistB(index).name];
  imgAinfo = imfinfo(fnameA);
  record = []; %variable to store each image data
  GUV = []; %variable to track each GUV
  disp(['File ', [num2str(index)],' open...']);
  
  for indexStack=(1+initSkip):1:size(imgAinfo,1)
    
  imgA = imread(fnameA, indexStack);
  imgB = imread(fnameB, indexStack);
  imgA = uint16(imgA); %casting added to avoid data type confusion in thresholding
  imgB = uint16(imgB);
  
  %% Binary mask created for high accuracy CHT detection 
  if autoBinaryMask == 1
    %maskThres = 65535 * graythresh(imgA, 'moments'); %auto moments thresholding
    
    [lehisto x] = imhist(imgA,11107); %triangle thresholding  
    maskThres = 65535 * triangle_th(lehisto, 11107);
    
  endif
  if autoBinaryMask == 0
     imgAbinary = imgA; %no ibnary processing
  else
     imgAbinary = imgA;
     imgAbinary(imgAbinary>maskThres) = 1;
     imgAbinary(imgAbinary != 1) = 0;
     imgAbinary = logical(imgAbinary);
  endif

  
  %%%%%%%%%%
    
  %CHT GUV pick up based on ChA
  figure(indexStack)
  [centers, radii, strengths] = imfindcircles (imgAbinary, [radius_min radius_max], 'sensitivity',sensitivity, 'edgethreshold',edge);
  %0.3 maximum image is more visible than full range images
  imshow(imgA, [min(min(imgA)) 0.5*max(max(imgA))])
  title(['ChA Stack image ',num2str(indexStack),' circle detection, Cyan-detected / Yellow-Filtered']);
  viscircles(centers, radii,'color','cyan');
     
  %background noise level analysis
  if autoBackground == 1
    backgroundIntA = mode(imgA(:));
    backgroundIntB = mode(imgB(:));
  elseif autoBackground == 2
    backgroundIntA = median(imgA(:));
    backgroundIntB = median(imgB(:));
  endif
  
  %Intensity analysis
  #[total_intensityA     number_of_pixelsA       average_intensityA
  # average_backgroundA   int-backgroundA      
  #total_intensityB     number_of_pixelsB       average_intensityB
  # average_backgroundB   int-backgroundB   
  # normbalizedB/A]
  intensity = zeros(length(radii),11);
  for i=1:1:size(intensity,1)
    r15 = radii(i) + 2 * thres; % little larger than radius
    for ix=round(centers(i,1)-r15):1: round(centers(i,1)+r15)
      for iy=round(centers(i,2)-r15):1: round(centers(i,2)+r15)
        if (ix<= 0) || (iy<=0) ||(ix>size(imgA,2))||(iy>size(imgB,1))
          continue
        endif
        dist2 = (ix - centers(i,1))^2 + (iy - centers(i,2))^2;
        if ((radii(i) - sqrt(dist2)) > 0)&&((radii(i) - sqrt(dist2)) < thres)
          intensity(i,1) += int32(imgA(iy, ix));    
          intensity(i,2) += 1;      
          intensity(i,6) += int32(imgB(iy, ix));    
          intensity(i,7) += 1;  
          
        endif
      endfor
    endfor
    intensity(i,3) = intensity(i,1) / intensity(i,2);
    intensity(i,4) = backgroundIntA; 
    intensity(i,5) = intensity(i,3) - backgroundIntA; 
    intensity(i,8) = intensity(i,6) / intensity(i,7);
    intensity(i,9) = backgroundIntB; 
    intensity(i,10) = intensity(i,8) - backgroundIntA; 
    %intensity(i,11) = intensity(i,10) / intensity(i,5); --removed as it is never used--
  endfor
  
disp(['Stack ', [num2str(indexStack)],' : total ',[num2str(length(radii))], ' circles were detected']);  
  
  
  %GUV filtering Strategies 1 (each image level) %%%%%%%%%%
  %1. Removing GUVs trapped in another GUV
  %% -- removed due to redundancy-- %%

  %2. MLV or intensity full image removal (or bottom stuck)
  if isMLVfilt == 1
    MLVCount = 0;
    IntGUVCount = 0;
    tempCenters = [];
    tempRadii = [];
    tempStrengths = [];
    tempIntensity = [];
    filteredCenters = [];
    filteredRadii = [];
    for i=1:1:size(intensity,1)
      r05 = radii(i) * 0.5;
      internalInt = 0;
      area = 0;
      for ix=round(centers(i,1)-r05):1: round(centers(i,1)+r05)
         for iy=round(centers(i,2)-r05):1: round(centers(i,2)+r05)
          if (ix<= 0) || (iy<=0) ||(ix>size(imgA,2))||(iy>size(imgB,1))
          continue
          endif
           internalInt += int32(imgA(iy, ix)); 
           area += 1;
         endfor
      endfor
      % outer rim intensity should be at least MLVthres x 100% of the inner intensity
      if intensity(i,5) < MLVthres * (internalInt/area - intensity(i,4))  
        %record{indexStack}(i-MLVCount,:) = []; %remove the image
        filteredCenters = [filteredCenters; centers(i,:)];
        filteredRadii = [filteredRadii; radii(i,:)];
        MLVCount += 1;
      else
        tempCenters =[tempCenters; centers(i,:)];
        tempRadii = [tempRadii; radii(i,:)];
        tempStrengths = [tempStrengths; strengths(i,:)];
        tempIntensity = [tempIntensity ; intensity(i,:)]; 
      endif
    endfor
    centers = tempCenters;
    radii = tempRadii;
    strengths = tempStrengths;
    intensity = tempIntensity;
    if length(filteredCenters) != 0
      viscircles(filteredCenters, filteredRadii, 'color','yellow','linestyle','--');
    endif
    disp([num2str(MLVCount),' GUVs removed for high intensity within the circle...']);
  endif
  
  record{indexStack} = [centers, radii, strengths, intensity];
  
  
  if saveCircles == 1
    print(indexStack, [PathOut '/' strrep(filelistA(index).name,'.tif',['_' num2str(indexStack) '.jpg'])]);  
  endif
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  endfor % stack loop
  
  %%%%%GUV stack tracking%%%%%%
  %%GUV info: [center x, center y, radii, strengths, intensity info..., stack#];
  
  for i=1:1:size(record{1},1)
   GUV{i} = [record{1}(i,:) 1]; % 1 at the end is stack number index information 
  endfor
        
  for indexStack=2:1:size(imgAinfo,1)
    for i=1:1:size(record{indexStack},1)
      ifFound = 0 ; %variable to decide if matching GUV was found 
      for j = 1:1:length(GUV)
        GUVchain = size(GUV{j},1); %length of tracking history for the GUV
        distGUV = sqrt((GUV{j}(GUVchain,1)-record{indexStack}(i,1))^2 + (GUV{j}(GUVchain,2)-record{indexStack}(i,2))^2);
        
        if ifFound == 1
          %nothing
        elseif distGUV < trackDist
          GUV{j} = [GUV{j}; record{indexStack}(i,:) indexStack];
          ifFound = 1;
        endif
        
      endfor
      if ifFound == 0 %start a new GUV coordinate if no match found
        GUV{ size(GUV,2)+1} = [record{indexStack}(i,:) indexStack]; 
      endif
    endfor
  endfor
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  disp([num2str(length(GUV)), ' GUV chains idenfied total before filtering...']);
  
  %GUV filtering Strategies 2 (stack level) %%%%%%%%%%
  % 1. Removing GUVs with too small number of z-Stacks
  if minNumStack > 0
    finalGUV = [];
    finalIndex = 1;
    filteredCount = 0;
    for indexGUV = 1:1:length(GUV)
      if size(GUV{indexGUV},1) >= minNumStack
       finalGUV{finalIndex} = GUV{indexGUV};
       finalIndex += 1;
      else
       filteredCount += 1;
      endif
    endfor
    GUV = finalGUV;
    disp([num2str(filteredCount),' GUVs removed for having less than ',num2str(minNumStack),' z-stacks...']);
  endif

  % 2. Removing GUVs that are stuck (too close in x,y space)
if (isCloseStackfilt == 1) && (length(GUV) > 1)
    finalGUV = [];
    finalIndex = 1;
    filteredCount = 0;
    for i = 1:1:length(GUV)
      ifFiltered = 0;
      for j = 1:1:length(GUV)
        if i == j
          continue;
        endif
        dist2 = (mean(GUV{i}(:,1)) - mean(GUV{j}(:,1)))^2 + (mean(GUV{i}(:,2)) - mean(GUV{j}(:,2)))^2;
        dist = sqrt(dist2);
        if dist < (max(GUV{i}(:,3)) +(max(GUV{j}(:,3))))
          filteredCount += 1;
          ifFiltered = 1;
        endif      
      endfor
        if ifFiltered == 0
         finalGUV{finalIndex} = GUV{i};
         finalIndex += 1; 
        endif
    endfor
  GUV = finalGUV;
  disp([num2str(filteredCount),' GUVs removed for too close proximity...']);
elseif isCloseStackfilt == 2

    finalGUV = [];
    finalIndex = 1;
    filteredCount = 0;
    for i = 1:1:length(GUV)
      ifFiltered = 0;
      for j = 1:1:length(GUV)
        if i == j
          continue;
        endif
        dist2 = (mean(GUV{i}(:,1)) - mean(GUV{j}(:,1)))^2 + (mean(GUV{i}(:,2)) - mean(GUV{j}(:,2)))^2;
        dist = sqrt(dist2);
        if (dist < (max(GUV{i}(:,3)) +(max(GUV{j}(:,3))))) && (mean(GUV{i}(:,4)) < mean(GUV{j}(:,4)))
          filteredCount += 1;
          ifFiltered = 1;
        endif      
      endfor
        if ifFiltered == 0
         finalGUV{finalIndex} = GUV{i};
         finalIndex += 1; 
        endif
    endfor
  GUV = finalGUV;
  disp([num2str(filteredCount),' GUVs removed for too close proximity...']);  
 
endif
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% Remove GUVs falsely identified due to lack of GUVs %%%%%%%%%%
  % When there is no GUVs to identify, CHT may pick up false one at the center, 
  % this part is to remove that x=(lengthX+1)/2, y=(lengthY+1)/2
  if minNumStack > 0
    finalGUV = [];
    finalIndex = 1;
    filteredCount = 0;
    for indexGUV = 1:1:length(GUV)
      if abs((mean(GUV{indexGUV}(:,2)) - (size(imgA,1) +1 )/2)) < 0.00000001 && abs((mean(GUV{indexGUV}(:,1)) - (size(imgA,2) +1 )/2)) <0.00000001    
       filteredCount += 1;
      else
       finalGUV{finalIndex} = GUV{indexGUV};
       finalIndex += 1;
      endif
    endfor
    GUV = finalGUV;
    disp([num2str(filteredCount),' GUVs removed for false indentification for no GUVs ']);
  endif
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%Save relevant information%%%%
  toSave = []; % matrix to save as txt
  for indexGUV = 1:1:length(GUV)
    avg_coordinate = [mean(GUV{indexGUV}(:,1)),mean(GUV{indexGUV}(:,2))];
    max_radii = max(GUV{indexGUV}(:,3));
    num_slices = size(GUV{indexGUV},1); 
    %intensity total A and B here to save 
    total_intensityA = sum(GUV{indexGUV}(:,5));
    number_of_pixelsA = sum(GUV{indexGUV}(:,6));
    average_intensityA = total_intensityA / number_of_pixelsA;
    backgrount_intensityA = sum (GUV{indexGUV}(:,8) .* GUV{indexGUV}(:,6)) /number_of_pixelsA;
    adjusted_intensityA = average_intensityA - backgrount_intensityA;
    
    total_intensityB = sum(GUV{indexGUV}(:,10));
    number_of_pixelsB = sum(GUV{indexGUV}(:,11));
    average_intensityB = total_intensityB / number_of_pixelsB;
    backgrount_intensityB = sum (GUV{indexGUV}(:,13) .* GUV{indexGUV}(:,11)) /number_of_pixelsB;
    adjusted_intensityB = average_intensityB - backgrount_intensityB;
    
    normalizedBA = adjusted_intensityB / adjusted_intensityA;

 
    GUV{indexGUV};
    
    toSave = [toSave; avg_coordinate max_radii num_slices total_intensityA number_of_pixelsA average_intensityA backgrount_intensityA adjusted_intensityA total_intensityB number_of_pixelsB average_intensityB backgrount_intensityB adjusted_intensityB normalizedBA];
  endfor
  
  OutName = strrep(filelistA(index).name,'.tif','.txt');
  dlmwrite([PathOut '\' OutName],toSave,'delimiter','\t')  
  if length(toSave) != 0
    disp(['Average intensity ChA: ', num2str(mean(toSave(:,7))), ' Average intensity ChB: ', num2str(mean(toSave(:,12)))]);
    disp(['Average adjusted intensity ChA: ', num2str(mean(toSave(:,9))), ' Average adjusted intensity ChB: ', num2str(mean(toSave(:,14)))]);
    disp(['Normalized intensity ChB/ChA (total total pixels): ', num2str(mean(toSave(:, 15)))]);
  endif
  
  disp([' '])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%Graphically show the final fitting result%%%%
  if length(toSave)>0
 
  imageToUse = round(indexStack/4);%use one image at the middle
  imgA = imread(fnameA, imageToUse);
  imgB = imread(fnameB, imageToUse);

  figure(indexStack + 1)
  imshow(imgA, [min(min(imgA)) 0.5*max(max(imgA))])
  title('ChA overlapped with final GUVs detected');
  centers = toSave(:,1:2); 
  radii = toSave(:,3);   
 
  if length(centers) != 0 
    viscircles(centers, radii, 'color','cyan')
  endif
  if saveCircles == 1
    print(indexStack + 1,[PathOut '/' strrep(filelistA(index).name,'.tif','_final_chA.jpg')]);  
  endif
  
  figure(indexStack + 2)
  imshow(imgB, [min(min(imgB)) 0.5*max(max(imgB))])
  title('ChB overlapped with final GUVs detected');
  if length(centers) != 0 
    viscircles(centers, radii, 'color','cyan')
  endif
  if saveCircles == 1
    print(indexStack + 2,[PathOut '/' strrep(filelistA(index).name,'.tif','_final_chB.jpg')]);  
  endif
  
     
  endif
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endfor





