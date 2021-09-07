%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUV phase state analysis using deep learning-based
% This program uses user defined variables (size of vesicles, sensitivity)
% to automatically pick up GUV, and analyze fluoresence intensity around it for
% phase state decision.
% Definition of macro-scale phase separated state:
% In a single image, at least 25% of the length from a detected circle should be
% in a distinct intensity or a phase region. Only when 50% or higher number of 
% image stakcs of a GUV is multi or binary phase separated it is called
% multi or binary phase separated. Otherwise just uniform. 
% In this neural network version, deep learning
% network trained on thousdands of images will be used to make an instant
% decision indpendent of the intensity profile. 
%
% Il-Hyung Lee (leei@montclair.edu), 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instruction: This Code was written for Matlab R2020b. Running it in other
% programs such as GNUOctave may need some modification. 
% Users can choose the target folder contianing the images, and 
% the saved neural network to use in .m format. triangle_th.m 
% should be accessible for this program. Images should
% be in Tiff format. Output files will be automatically generated. 
% Output images show initial circle detection for each z-section and then
% state decisions. Stack based final output + each intermediate z-stack image
% calculation are also saved as well for any intensity analysis. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pkg load image; %legacy command for GNUOctave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output txt format 
%% Each raw is for one GUV detected
%% from left to right:
%% 1. Average coordinate x of the GUV
%% 2. Average coordinate y of the GUV
%% 3. Maximum radius of the GUV (Greatest raadius of multiple Z-images)
%% 4. Total number of Z-images of the GUV
%% 5. Average per pixel intensity of the GUV
%% 6. minimum per pixel intensity
%% 7. maximum per pixel intensity 
%% 8. Phase decision (0-uniform, 1- binary pahse, 2- more than binary phase by definition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radius_min = 20; ; %smallest radius of GUV to detect
radius_max = 120; ; %largest radius of GUV to detect
sensitivity = 0.90; %CHT fit sensitivity (0-1, 1 means free pass)
edge = 0.1; %CHT edge definition sensitivity (0-1, 1 means trong edge needed)
thres = 10;%4; %number of pixels to count is within the fluorescence
thres2 = thres^2; 
num_seg = 72; %umber of segments to devide polar coordinate of GUVs should
% be a number can can be devided by 2 
trackDist = 5; %pixel distance to consider the same vesicles

%%%%%%%%%%
initSkip = 0; %%number of images to skip (in case initial few are not good)


autoBackground = 1; %automatic estimation of background intensity by most probable 
                    %vlaue 1 to use mode estimation, value 2 for median estimation
backgroundIntA = 0;  %manual input of the backgroudn intensity in case not using auto estimation
%backgroundIntB = 0;

saveCircles = 1; %1 - on to save all circle figures as jpeg 
saveEachIntensity = 1; %1 - on to save all circle intensity history (including GUV filtered ones)

phasDecision = 0.1; %variable to make decision about discontinuity around fluorescence contour
                    %It mens Lowintensity <= Highintensity x 200% x phasEcision  
                    %determine this based on manual examination of per pixel intentisy of two
                    %different phases
decisionP = 0.4; %Vesicles with decisionPx100% of phase separation will be called phase separated. 

                    
autoBinaryMask = 1; % 1 to use automatic binary mask thesholding 2 to use manual mask value 0 to use raw images
maskThres = 0.5; % if not using automask, any intensity greater than this avlue will be 1 and smaller will be 0

minNumStack = 5; % Minimum number of z-stacks needed to count as GUV (0 for not using it)
ifFilterEdgeVes = 1; %1 means vesicles partially detected at the edge will be excluded for analysis. 
                     % This should be on unless necessary. 
isNRfilt = 1; % 1 for using neural decision to filter invalid GUV images

isCloseStackfilt = 1; % 1 for filtering GUV Stacks that are too close
                      % 1 to remove both GUVs stuck 2 to reject only one with less fitting fitness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PathName=uigetdir('Choose file(stack) containing directory (single color)');
PathName = './Sample/';
filelistA = dir([PathName, '/*.tif']);
PathOut = [PathName '/Analysis'];
if exist([PathName '/Analysis']) ~= 7
  mkdir(PathOut);
end

%%specify the trained neural network to use%%
%[neuralFile neuralPath] = uigetfile('*.mat','Choose the neural network file to load.');
neuralPath = './';
neuralFile = 'nrFilter4C.mat';
load([neuralPath, neuralFile]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=1:1:length(filelistA)
 
  fnameA = [PathName '\' filelistA(index).name];
  imgAinfo = imfinfo(fnameA);
  record = []; %variable to store each image data
  GUV = []; %variable to track each GUV
  
  disp(['File ', [num2str(index)],' open...']);
  
  for indexStack=(1+ initSkip):1:size(imgAinfo,1)
      
  imgA = imread(fnameA, indexStack);
  imgA = uint16(imgA); %casting added to avoid data type confusion in thresholding
  
  %% Binary mask created for high accuracy CHT detection 
  if autoBinaryMask == 1
    %maskThres = 65535 * graythresh(imgA, 'moments'); %auto moments thresholding
    
    [lehisto x] = imhist(imgA,11107); %triangle thresholding  
    maskThres = 65535 * triangle_th(lehisto, 11107);
    
  end
  if autoBinaryMask == 0
     imgAbinary = imgA; %no ibnary processing
  else
     imgAbinary = imgA;
     imgAbinary(imgAbinary>maskThres) = 1;
     imgAbinary(imgAbinary ~= 1) = 0;
     imgAbinary = logical(imgAbinary);
  end
  
  %CHT GUV pick up
  figure(indexStack)
  [centers, radii, strengths] = imfindcircles (imgAbinary, [radius_min radius_max], 'sensitivity',sensitivity, 'edgethreshold',edge);
  imshow(imgA, [min(min(imgA)) 0.5 * max(max(imgA))])
  title(['ChA Stack image ',num2str(indexStack),' circle detection Cyan-valid Red-edge Green-C2 Blue-C3 Magenta-C4 Yellow C5']);
  viscircles(centers, radii,'color','cyan');
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %background noise level analysis
  if autoBackground == 1
    backgroundIntA = mode(imgA(:));
    %backgroundIntB = mode(imgB(:));
  elseif autoBackground == 2
    backgroundIntA = median(imgA(:));
    %backgroundIntB = median(imgB(:));
  end
  
  %segmented Intensity analysis for phase decision
  %[total_intensity total_pixels avg_intensity backgroundIntA Adjusted_intensity
  % adjIntensity_Seg1 2 3 4 5 ... 18 phase_decision Highintadj Lowintadj]
  intensity = zeros(length(radii),3+2 + num_seg + 3);
 
  
  for i=1:1:size(intensity,1)
    seg_pixels = zeros(1, num_seg); %to store number of pixels involves in calc
    r15 = radii(i) + 2 * thres; % little larger than radius
    rsqrd = radii(i)^2; %radius squared
    for ix=round(centers(i,1)-r15):1: round(centers(i,1)+r15)
      for iy=round(centers(i,2)-r15):1: round(centers(i,2)+r15)
        if (ix<= 0) || (iy<=0) ||(ix>size(imgA,2))||(iy>size(imgA,1))
         continue
        end
        ix_coord = ix - centers(i,1); %coordinate relative to the centers
        iy_coord = iy - centers(i,2);
        dist2 = (ix_coord)^2 + (iy_coord)^2;
        dist = sqrt(dist2);
        if ((radii(i) - dist) > 0)&&((radii(i) - dist) < thres)
        %%if (abs(radii(i) - dist) < thres)
          intensity(i,1) = intensity(i,1) + int32(imgA(iy, ix));    
          intensity(i,2) = intensity(i,2) + 1;      
          %%%segement analysis%%%
          if iy_coord >= 0  % coordinate y>0
            seg = 1 + floor(num_seg/2 * (acos(ix_coord/dist)/pi));
            intensity(i, 3+2+seg)=intensity(i, 3+2+seg) + int32(imgA(iy, ix));
            seg_pixels(seg) = seg_pixels(seg) + 1;            
          else % coordinate y<0
            seg = num_seg - floor(num_seg/2 * (acos(ix_coord/dist)/pi));
            intensity(i, 3+2+seg)= intensity(i, 3+2+seg) + int32(imgA(iy, ix));
            seg_pixels(seg) = seg_pixels(seg) + 1;    
          end
          
        end
      end
    end
    
    for seg = 1:1:num_seg
      if seg_pixels(seg) == 0
        intensity(i, 3+2+seg) = 0;
      else
        intensity(i, 3+2+seg) = intensity(i, 3+2+seg) / seg_pixels(seg);
        intensity(i, 3+2+seg) = intensity(i, 3+2+seg) - backgroundIntA;
      end
    end 
    intensity(i,3) = intensity(i,1) / intensity(i,2);
    intensity(i,4) = backgroundIntA;
    intensity(i,5) = intensity(i,3) - backgroundIntA;
  end
  disp(['Stack ', [num2str(indexStack)],' : total ',[num2str(length(radii))], ' circles were detected']);
  
  %%%%%%Single image based filtering%%%%%  
  %GUV filtering Strategies 1 (each image level) %%%%%%%%%%
  %0. Remove GUVs on the edges (incomplete circles)
  if ifFilterEdgeVes == 1
    tempCenters = [];
    tempRadii = [];
    tempStrengths = [];
    filteredCenters = [];
    filteredRadii = [];
    tempIntensity = [];
    
    for i=1:1:length(radii)
      if (centers(i, 1) + radii(i)) > size(imgA,2) || (centers(i, 1) - radii(i)) < 0 || (centers(i, 2) + radii(i)) > size(imgA,1) || (centers(i, 2) - radii(i)) < 0
      %bypass the vesicle
      filteredCenters = [filteredCenters; centers(i,:)];
      filteredRadii = [filteredRadii; radii(i,:)];
      else
      tempCenters =[tempCenters; centers(i,:)];
      tempRadii = [tempRadii; radii(i,:)];
      tempStrengths = [tempStrengths; strengths(i,:)];
      tempIntensity = [tempIntensity ; intensity(i,:)]; 
      
      end
    end
    centers = tempCenters;
    radii = tempRadii;
    strengths = tempStrengths;
    intensity = tempIntensity;
    
    
    if length(filteredCenters) ~= 0
      viscircles(filteredCenters, filteredRadii, 'color','red','linestyle','--');
    end
  
  disp([num2str(length(filteredRadii)),' Circles removed for being at the edges...']);  
  end
  
  %1. neural decision filter to remove invalid GUVs
  if size(intensity,1) > 1
  if isNRfilt == 1
    IntGUVCount = 0;
    tempCenters = [];
    tempRadii = [];
    tempStrengths = [];
    tempIntensity = [];
    
    filteredCenters2 = [];
    filteredRadii2 = [];
    filteredCenters3 = [];
    filteredRadii3 = [];
    filteredCenters4 = [];
    filteredRadii4 = [];
    filteredCenters5 = [];
    filteredRadii5 = [];
    
    for i=1:1:size(intensity,1)
  
        nrMargin = 1.2; % to multiply to make sure image contains the circle.
        nrX1 = round(centers(i,1) - radii(i)*nrMargin); 
        if nrX1 <= 0
            nrX1 = 1;
        end
        nrX2 = round(centers(i,1) + radii(i)*nrMargin); 
        if nrX2 > size(imgA,2)
            nrX2 = size(imgA,2);
        end
            nrY1 = round(centers(i,2) - radii(i)*nrMargin); 
        if nrY1 <= 0
            nrY1 = 1;
        end
            nrY2 = round(centers(i,2) + radii(i)*nrMargin); 
        if nrY2 > size(imgA,1)
            nrY2 = size(imgA,1);
        end
        
        nrImage = imgA(nrY1:nrY2,nrX1:nrX2);
        nrImage = imresize(nrImage,[50, 50]);
        
        [decision, probNR] = classify(net,nrImage);
      
          %This part can take from 2 to 5 classes
          if decision == 'C2'
            IntGUVCount = IntGUVCount + 1; 
            filteredCenters2 = [filteredCenters2; centers(i,:)];
            filteredRadii2 = [filteredRadii2; radii(i,:)];
          
          elseif decision == 'C3'
            IntGUVCount = IntGUVCount + 1;  
            filteredCenters3 = [filteredCenters3; centers(i,:)];
            filteredRadii3 = [filteredRadii3; radii(i,:)];
          
          elseif decision == 'C4'
            IntGUVCount = IntGUVCount + 1;
            filteredCenters4 = [filteredCenters4; centers(i,:)];
            filteredRadii4 = [filteredRadii4; radii(i,:)];
          elseif decision =='C5'
            IntGUVCount = IntGUVCount + 1;
            filteredCenters5 = [filteredCenters5; centers(i,:)];
            filteredRadii5 = [filteredRadii5; radii(i,:)];  
                    
          else
            tempCenters =[tempCenters; centers(i,:)];
            tempRadii = [tempRadii; radii(i,:)];
            tempStrengths = [tempStrengths; strengths(i,:)];
            tempIntensity = [tempIntensity ; intensity(i,:)]; 
            
          end           
      
    end
    centers = tempCenters;
    radii = tempRadii;
    strengths = tempStrengths;
    intensity = tempIntensity;
    
    
      viscircles(filteredCenters2, filteredRadii2, 'color','green','linestyle','--');
    
      viscircles(filteredCenters3, filteredRadii3, 'color','blue','linestyle','--');
   
      viscircles(filteredCenters4, filteredRadii4, 'color','magenta','linestyle','--');
    
      viscircles(filteredCenters5, filteredRadii5, 'color','yellow','linestyle','--');  
    
  disp([num2str(IntGUVCount),' Circles removed for being invalid...']); 
  
  end
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
  %%%%%%%%phase decision%%%%%%%%%%%%%%%%%%%%%%%
 
  
  if size(intensity,1) ~= 0
  for i=1:1:size(intensity,1)  
    %intMax = max(intensity(i,6:5+num_seg));
    %intMin = min(intensity(i,6:5+num_seg));
    sortIntensity = sort(intensity(i,6:5+num_seg));
    intMax = sortIntensity(round(0.8*num_seg));
    %intMax is an estimated intensity for a high intensity domain 
    %it is estimated as a pixel intensity of 20% from the highest
    intMin = sortIntensity(round(0.2*num_seg));
    %intMin is an estimated intensity for a low intensity domain 
    %it is estimated as a pixel intensity of 20% from the lowest
    
    
    intMean = (intMax + intMin)/2;
    intMeanHigh = intMean + intMean * phasDecision;
    intMeanLow = intMean - intMean * phasDecision;
    temp = 0; %temporary variable to judge current intensity level
    phase = 0; %total number of phases
    intHigh = [];
    intLow = [];
    if intensity(i,5+num_seg) > intMeanHigh
      temp = 1;
    elseif intensity(i,5+num_seg) < intMeanLow
      temp = -1;
    else
      temp = 0;
    end
    
    for j=1:1:num_seg
      if temp == 1
        if intensity(i,5+j) > intMeanHigh
          %nothing
          intHigh = [intHigh intensity(i,5+j)];
        elseif intensity(i,5+j) < intMeanLow
          temp = -1;
          phase = phase + 1;
          intLow = [intLow intensity(i,5+j)];
        else %intensity medium
          %nothing
        end
          
      elseif temp == -1
        if intensity(i,5+j) > intMeanHigh
          temp = 1;
          phase = phase + 1;
          intHigh = [intHigh intensity(i,5+j)];
        elseif intensity(i,5+j) < intMeanLow
          %nothing
          intLow = [intLow intensity(i,5+j)];
        else %intensity medium
          %nothing
        end    
      else %temp == 0 
        if intensity(i,5+j) > intMeanHigh
          temp = 1;
         elseif intensity(i,5+j) < intMeanLow
          temp = -1;
        else %intensity medium
          %nothing
        end    
      end
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    %% final phase decision
    if (length(intHigh)/num_seg > 0.25 ) || (length(intLow)/num_seg > 0.25 )
      phase = phase;
    else
      phase = 0;
    end
    %%
    
    intensity(i,5+num_seg+1) = phase/2;
    if length(intHigh)>0
      intensity(i,5+num_seg+2) = median(intHigh);
    else
      intensity(i,5+num_seg+2) = 0;
    end
    
    if length(intLow)>0
      intensity(i,5+num_seg+3) = median(intLow);
     else
      intensity(i,5+num_seg+3) = 0;  
    end
    
  end
  end   
  
  

    record{indexStack} = [centers, radii, strengths, intensity]; %%NR
  
  if saveCircles == 1
    print(indexStack, [PathOut '/' strrep(filelistA(index).name,'.tif',['_' num2str(indexStack) '.jpg'])], '-djpeg');  
  end
  if saveEachIntensity == 1
    OutName = strrep(['Int_' filelistA(index).name],'.tif',['_' num2str(indexStack) '.txt']);
    dlmwrite([PathOut '\' OutName],record{indexStack},'delimiter','\t');  
  end
  
  end % stack loop
  
  %%%%%GUV stack tracking%%%%%%
  
  for i=1:1:size(record{1},1)
   GUV{i} = record{1}(i,:);
  end
        
  for indexStack=2:1:size(imgAinfo,1)
    for i=1:1:size(record{indexStack},1)
      ifFound = 0 ; %variable to decide if matching GUV was found 
      for j = 1:1:length(GUV)
        GUVchain = size(GUV{j},1); %length of tracking history for the GUV
        distGUV = sqrt((GUV{j}(GUVchain,1)-record{indexStack}(i,1))^2 + (GUV{j}(GUVchain,2)-record{indexStack}(i,2))^2);
        
        if ifFound == 1
          %nothing
        elseif distGUV < trackDist
          GUV{j} = [GUV{j}; record{indexStack}(i,:)];
          ifFound = 1;
        end
        
      end
      if ifFound == 0 %start a new GUV coordinate if no match found
        GUV{ size(GUV,2)+1} = record{indexStack}(i,:); 
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  disp([num2str(length(GUV)), ' GUV chains idenfied total before filtering...']);
  
  %%%%%GUV filtering%%%%%%%%%%%
  %GUV filtering Strategies 2 (stack level) %%%%%%%%%%
  % 1. Removing GUVs with too small number of z-Stacks
  if minNumStack > 0
    finalGUV = [];
    finalIndex = 1;
    filteredCount = 0;
    for indexGUV = 1:1:length(GUV)
      if size(GUV{indexGUV},1) >= minNumStack
       finalGUV{finalIndex} = GUV{indexGUV};
       finalIndex = finalIndex + 1;
      else
       filteredCount = filteredCount + 1;
      end
    end
    GUV = finalGUV;
    disp([num2str(filteredCount),' GUVs removed for having less than ',num2str(minNumStack),' z-stacks...']);
  end
  
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
        end
        dist2 = (mean(GUV{i}(:,1)) - mean(GUV{j}(:,1)))^2 + (mean(GUV{i}(:,2)) - mean(GUV{j}(:,2)))^2;
        dist = sqrt(dist2);
        if dist < (max(GUV{i}(:,3)) +(max(GUV{j}(:,3))))
          filteredCount = filteredCount + 1;
          ifFiltered = 1;
        end      
      end
        if ifFiltered == 0
         finalGUV{finalIndex} = GUV{i};
         finalIndex = finalIndex + 1; 
        end
    end
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
        end
        dist2 = (mean(GUV{i}(:,1)) - mean(GUV{j}(:,1)))^2 + (mean(GUV{i}(:,2)) - mean(GUV{j}(:,2)))^2;
        dist = sqrt(dist2);
        if (dist < (max(GUV{i}(:,3)) +(max(GUV{j}(:,3))))) && (mean(GUV{i}(:,4)) < mean(GUV{j}(:,4)))
          filteredCount = filteredCount + 1;
          ifFiltered = 1;
        end      
      end
        if ifFiltered == 0
         finalGUV{finalIndex} = GUV{i};
         finalIndex = finalIndex + 1; 
        end
    end
  GUV = finalGUV;
  disp([num2str(filteredCount),' GUVs removed for too close proximity...']);  
 
end
  
  
  %% 3. Remove GUVs falsely identified due to lack of GUVs %%%%%%%%%%
  % When there is no GUVs to identify, CHT will pick up false one at the center, 
  % this part is to remove that x=(lengthX+1)/2, y=(lengthY+1)/2
  if minNumStack > 0
    finalGUV = [];
    finalIndex = 1;
    filteredCount = 0;
    for indexGUV = 1:1:length(GUV)
      if abs((mean(GUV{indexGUV}(:,2)) - (size(imgA,1) +1 )/2)) < 0.00000001 && abs((mean(GUV{indexGUV}(:,1)) - (size(imgA,2) +1 )/2)) <0.00000001    
       filteredCount = filteredCount + 1;
      else
       finalGUV{finalIndex} = GUV{indexGUV};
       finalIndex = finalIndex + 1;
      end
    end
    GUV = finalGUV;
    disp([num2str(filteredCount),' GUVs removed for false indentification for no GUVs ']);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%Save relevant information%%%%
  %file number, GUV number, avg_coordinate, max_radii, #of slices, 
  % avg fluo intensity, max fluo, min fluo, phase decision, phase%, phase info?
  toSave = []; % matrix to save as txt
  countUniform = 0;
  countPhase = 0;
  countMulti = 0;

  
  for indexGUV = 1:1:length(GUV)
    avg_coordinate = [mean(GUV{indexGUV}(:,1)),mean(GUV{indexGUV}(:,2))];
    max_radii = max(GUV{indexGUV}(:,3));
    num_slices = size(GUV{indexGUV},1); 
    avg_intensity = mean(GUV{indexGUV}(:,7));
    max_intensity = max(max(GUV{indexGUV}(:,8:7+num_seg)));
    min_intensity = min(min(GUV{indexGUV}(:,8:7+num_seg)));
    
    % if more than decisionP % is phase separated, it is defined as phase separated
    sortGUV = sort(GUV{indexGUV}(:,2+8+num_seg));
    if sortGUV(round((1-decisionP)*length(sortGUV))) >= 2
      phaseGUV = 2; %multi phase
      countMulti = countMulti + 1;
    elseif sortGUV(round((1-decisionP)*length(sortGUV))) >=1
      phaseGUV = 1; %two phases
      countPhase = countPhase + 1;
    else  
      phaseGUV = 0; %uniform
      countUniform = countUniform + 1;
    end
    

    toSave = [toSave; avg_coordinate max_radii num_slices avg_intensity max_intensity min_intensity phaseGUV];
  end
  
  disp(['Total ', num2str(length(GUV)),' GUVs. Uniform: ', num2str(100*countUniform/length(GUV)),'%, Phase: ', num2str(100*countPhase/length(GUV)),'%, Multi: ', num2str(100*countMulti/length(GUV)), '%']);
  disp(' ');
  
  OutName = strrep(filelistA(index).name,'.tif','.txt');
  dlmwrite([PathOut '\' OutName],toSave,'delimiter','\t') ; 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%Graphically show the final fitting result%%%%
  if length(toSave)>0
 
  if size(imgAinfo,1) >4
   imageToUse = 5;%use one image at the middle
   imgA = imread(fnameA, imageToUse);
   
  else
   imgA = imread(fnameA,1);
  end

  figure(size(imgAinfo,1) + 1);
  imshow(imgA, [min(min(imgA)) 0.5*max(max(imgA))]);
  title('ChA overlapped with final GUVs detected. Cyan-uniform, Magenta-binary, Yellow-multi');
 
  sizeLength = size(toSave,2);
  
  for i = 1:1:size(toSave,1)
    if toSave(i,sizeLength) == 0
      viscircles(toSave(i,1:2), toSave(i,3), 'color','cyan','linestyle','-');
      
    elseif toSave(i,sizeLength) == 1
      viscircles(toSave(i,1:2), toSave(i,3), 'color','magenta','linestyle','--');
      
    else
      viscircles(toSave(i,1:2), toSave(i,3), 'color','yellow','linestyle',':');
     
    end
    
  end
  
  if saveCircles == 1
    print(size(imgAinfo,1) + 1,[PathOut '/' strrep(filelistA(index).name,'.tif','_final_chA.jpg')], '-djpeg');  
  end
 
  
  %%%%%%%%%%%%%%%%%%%%
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end




