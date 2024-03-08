
%%                                             IMAGE SEQUENCE TO MOVIE 2.0

%                                                Fabio Crameri, 29.09.2016
clear

%% INPUT
SAVE.FigureName         = 'test_fslip_eta*';                    %figure name - use * to account for sequential numbers
SAVE.Directory          = '/work/plotTest/+im/test_fslip/';     %directory of the figure sequence
SAVE.MovieName          = '+test';                              %movie name without format ending

SAVE.avi                = logical(1);
SAVE.mj2                = logical(0);
SAVE.mp4                = logical(0);
SAVE.m4v                = logical(0);

SAVE.FrameRate        	= 10;                                   %30 (default) | frames per second
SAVE.Quality          	= 75;                                   %75 (default) | integer in the range [0,100]

showImage               = logical(1);
SAVE.overwrite          = logical(0);
SAVE.blackBackground    = logical(0);

%--------------------------------------------------------------------------
% do not make changes to the code below this line

%% CHECK CODE STATUS
fAIo.app='I2M';
fAIo.appVersion=2.0;
fAIo.appDate='29-Sep-2016';
[fAIo] = f_AIoI2M(fAIo,'StartUp',SAVE);

%% MOVE TO DIRECTORY
cd(SAVE.Directory)
figureList = dir(SAVE.FigureName);
if size(figureList,1)==0
    disp(['No figure file "',SAVE.FigureName,'" found in ',SAVE.Directory])
else
    numFiguresFound = size(figureList,1);
    disp([num2str(numFiguresFound),' figure files found.'])
end

%% CREATE FRAMES
if showImage
   figure(1),clf
end
for ifig=1:numFiguresFound
    im = imread(figureList(ifig).name); %read image
    SAVE.MovieFrames(ifig) = im2frame(im);
    if showImage
       image(im) 
       axis equal
       axis off
    end
end

%% SAVE MOVIE
[fAIo] = f_AIoI2M(fAIo,'SavingMovie',SAVE);



