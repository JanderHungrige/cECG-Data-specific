% Reading in the annotation file for the cECG dataset
clc
clear


epoch=30; %[30] [1] 30 seconds or 1 secondepoch annotations
deleteOldMatfiles=1;

Annotator='B3A';
% Annotator='Va13ntina';

Neonate=[4,5,6,7,9,10,11,12,13];

for p=1:length(Neonate)
    pat=Neonate(p);
%     annotationfolder= ['E:\cECG_study\Annotations\' Annotator '\participant' num2str(pat)];
%     annotationfolder=['\\code1\storage\2012-0194_neonatal_data\cECG study\Annotations\Data and annotations used by Jan (Bea)\Annotation\participant' num2str(pat)];
    annotationfolder=['C:\Users\310122653\Documents\PhD\Article_3_(cECG)\Raw Data\Annotation\participant'  num2str(pat)];
    annotationfiles=dir([annotationfolder '\*' num2str(epoch) '.csv ' ]);
    savefolder=[annotationfolder '\'];
    
    if deleteOldMatfiles==1
        Matfiles=dir([annotationfolder '\*.mat']);
        for k=1:length(Matfiles)
            delete([annotationfolder '\' Matfiles(k).name])
        end
    end
    %% ***** READIN FILE *****
    % Read in header
    fileID=fopen([annotationfiles.folder '\' annotationfiles.name],'r');
    tline = fgetl(fileID);
    header(1,:) = regexp(tline, '\,', 'split');%  Split header

    %Readin annotations (with first line header)
    for ctr=2:9
        if ischar(tline)    
              tline = fgetl(fileID);         
              header(ctr,:) = regexp(tline, '\,', 'split'); 
        else
              break;     
        end
    end

    %  Parse and read rest of file
    ctr = 10;
    while(~feof(fileID))
        if ischar(tline)    
              ctr = ctr + 1;
              tline = fgetl(fileID);         
              annot(ctr-10,:) = regexp(tline, '\,', 'split'); 
        else
              break;     
        end
    end

    fclose(fileID);
    clearvars ctr same tline 

    % separate annotations per session
    annot_no_header=annot(2:end,:); % otherwise confusion with index later
    Sessionname=annot{2,1}; k=1;
    Sessionnames{k}=annot{2,1}; k=2; %collecting the first session name
    % collecting the rest of the session names
    for i=2:length(annot_no_header)
        same=strcmp(annot_no_header{i,1},Sessionname); % compare if the name is still the same
        if same==0 % if not...
            endsessions(k-1,1)=i-1; %write index in file... -1 because the cindex marks the first of new name not the end
            Sessionname=annot_no_header{i,1}; % and change Sessionname to the new filename
            Sessionnames{k}=Sessionname; %collecting all sessionames
            k=k+1;
        end
    end

    if exist('endsessions', 'var')==0 % if only one session exist, no endsessions created
        endsessions=length(annot_no_header);
    end

    %% ***** Creating sleep state array *****


    Annotations = zeros(length(annot_no_header),1);
    state ={...
        '		ActiveSleep';...
        '		QuietSleep';...
        '		Wake';...
        '		CareTaking';...
        '		UnknownBedState';...
        '		Transition'...
        };
    statevalue=[1;2;3;4;5;6];%12;21;13;31;23;32]; 6 = transition

    Position = zeros(length(annot_no_header),1);
    pos={'		Unknown';...
         '		Back'};
    posvalue=[8,9];

    %Transfer annotations from the excel file to an annotation array
    for i=1:length(state)
        Annotations=transferannot(annot_no_header,state{i},statevalue(i),Annotations);
%         Annotations=addTransition(annot_no_header,i,statevalue,Annotations); %change 6 into 61 62 63...         
    end
    
    %Transfer the positions annotation from the excel file to an annotation array
    for i=1:length(pos)
        Position=posAnnot(annot_no_header,pos{i},posvalue(i),Position);
    end

    %adapting name
    for i=1:length(Sessionnames)
        s = strsplit(Sessionnames{1,i},'.'); 
        Sessionnames{1,i}=s{2};
        Sessionpatient{1,i}=s{1};
%         Sessionnames{1,i}=regexprep(Sessionnames{1,i},{'\.','avi'},{''}); % delete .avi from name  
    end

    %Create annotations for each session of the patient
    for i=1:length(Sessionnames)
        if i==1 % For the first session            
            eval(['Annotations_' Sessionnames{1,i} '_' Sessionpatient{1,i} '= Annotations(1:endsessions(i),1);' ])
            eval (['Position_' Sessionnames{1,i} '_' Sessionpatient{1,i} '= Position(1:endsessions(i),1);' ])
        elseif i~=length(Sessionnames) % For any session inbetween 
            eval(['Annotations_' Sessionnames{1,i} '_' Sessionpatient{1,i} '= Annotations((endsessions(i-1)+1):endsessions(i),1);'])        
            eval (['Position_' Sessionnames{1,i} '_' Sessionpatient{1,i} '= Position((endsessions(i-1)+1):endsessions(i),1);'])            
        elseif i==length(Sessionnames) % For the last session
            eval(['Annotations_' Sessionnames{1,i} '_' Sessionpatient{1,i} '= Annotations((endsessions(i-1)+1):end,1);' ])                
            eval (['Position_' Sessionnames{1,i} '_' Sessionpatient{1,i} '= Position((endsessions(i-1)+1):end,1);'])            
            
        end
        
        eval(['save( "'  savefolder 'Annotations_' Sessionnames{1,i} '_' Sessionpatient{1,i} '","Annotations_' Sessionnames{1,i} '_' Sessionpatient{1,i}  '","Position_' Sessionnames{1,i} '_' Sessionpatient{1,i} '")'])      
    end
    clearvars Sessionnames Sessionpatient annot
end

%% Nested functions

% Transfering the annotations into one array with numbers for states
function Annotations=transferannot(annot,state,statevalue,Annotations)
    for i=1:length(annot)
        same=strcmp(annot{i,10},state);
        sameIS=strcmp(annot{i,11},state);
      if same==1
          Annotations(i,1)=statevalue;
      end
      if sameIS==1
          Annotations(i,1)=statevalue;
      end      
    end
end

%Changing Transition from AS to QS in 12 and QS to AS in 21
function Annotations=addTransition(annot,state,statevalue,Annotations)

% First determine just if transition. Later we need to deterine from
% which-to-which state.
    for i=1:length(annot)
        same=strcmp(annot{i,11},'		Transition');
      if same==1 && state==1 %AS
          Annotations(i,1)=statevalue(61,1);
      elseif  same==1 && state==2 %QS
          Annotations(i,1)=statevalue(62,1);      
      elseif  same==1 && state==3 %Wake
          Annotations(i,1)=statevalue(63,1);
      elseif  same==1 && state==4 %CareTaking
          Annotations(i,1)=statevalue(64,1);       
      elseif  same==1 && state==5 %Unknown Bedstate
          Annotations(i,1)=statevalue(65,1);              
      end
    end
end
% determine in which position the baby is
function Position=posAnnot(annot,pos,posvalue,Position)
    for i=1:length(annot)
        same=strcmp(annot{i,12},pos);
      if same==1
          Position(i,1)=posvalue;
      end        
    end
end


