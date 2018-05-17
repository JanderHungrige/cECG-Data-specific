%Splitting the Video and data in max 1 hour length
clc
clear
addpath('C:\Users\310122653\Documents\PhD\cECG Data\Quick Annotator')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Quick Annotator')
drive='e';
participant=3;
path=([drive ':\cECG_study\For_Quick_Annotator\participant' num2str(participant) ]);
camPath=([path '\camera\']);
addpath(camPath)
cd(camPath)
folder=dir('2012*');

%Loading ffmpeg into matlab

% path1 = getenv('PATH');
% path1 = [path1 'C:\ffmpeg\bin\'];
% setenv('PATH', path1);
% !echo $PATH 
% setenv PATH "${PATH}:C:\ffmpeg\bin\"
% getenv('PATH')

% VIDEO
numFiles=size(folder,1);
for i=1:numFiles
    cd([camPath folder(i,1).name])
    names=dir('*.avi');
    numVideos=size(names,1);
    for j=1:numVideos 
        info = mmfileinfo(names(j,1).name); %Getting length of Video
        HoursOfMovie=info.Duration/3600; % Get how many parts should be created
        if info.Duration >= 3600 %longer than 1 hour
            split=3600;
            if HoursOfMovie-floor(HoursOfMovie)<0.25 %If the duration of last video would be shorter than 15min use differnt split
                HoursOfMovie=info.Duration/2700; %split on 45min instead of 1 hour
                split=2700; %indication for name change to take 45min(2700s) instead of 1 hour (3600s)
            end           
            %reading Video
            vid=VideoReader([camPath folder(i,1).name '\' info.Filename]);
           
            for k=0:ceil(HoursOfMovie)-1 
                v.CurrentTime = 2.5;
%                 Film=vid.read([camPath folder(i,1).name '\' vid.Name],[k*split*vid.FrameRate k+split*vid.FrameRate]);
                %writing Video
                if split == 3600
                    add=k*3600;            
                elseif split ==2700
                    add=k*2700;   
                end
                s = strsplit(info.Filename,'.');                 
                name=[s{1,1} '.' num2str(str2num(s{1,2})+add)];  
                cd([path '\camera\' folder(i,1).name '\']) %go to videofolder
                system(['C:\ffmpeg\bin\ffmpeg.exe -i ' vid.Name ' -ss ' k*add ' -t ' k+1*add ' name.avi']);
                
%                 newvideo=VideoWriter([name '.avi']);                             
%                 open(newvideo);
%                 writeVideo(newvideo,Film)
           end

            
        end % if longer than 1 hour
    end % for several videos
end
