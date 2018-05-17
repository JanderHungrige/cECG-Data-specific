% Inter rater reliability

clear
clc


Annotators={'B3A';'Va13ntina'};

Neonate=[4,5,6,7,9,10,11,12,13];

        savefolder=['E:\cECG_study\Annotations\Interratervariability\'];

        
for p=1:length(Neonate)
pat=Neonate(p);

    for i=1:length(Annotators)
        annotationfolder= ['E:\cECG_study\Annotations\' Annotators{i} '\participant' num2str(pat)];
        annotationfiles=dir([annotationfolder '\*.mat ' ]);
            %adapting name
        for p=1:length(annotationfiles)
            names{p,1}=regexprep(annotationfiles(p,1).name,{'\.','mat'},{''}); % delete .avi from name  
        end
        
        
        for m=1:length(annotationfiles)
            tmp=load([annotationfolder '\' annotationfiles(m).name]);       
            eval([names{m,1} '_' Annotators{i} '= tmp.' names{m,1} ';'])
            clearvars tmp
        end
    end
%% Kappa        
% Syntax: 	kappa(X,W,ALPHA)
%      
%     Inputs:
%           X - square data matrix
%           W - Weight (0 = unweighted; 1 = linear weighted; 2 = quadratic
%           weighted; -1 = display all. Default=0)
%           ALPHA - default=0.05.        
for n=1:length(annotationfiles)
    eval(['M1(1,:)=' names{m,1} '_' Annotators{1} ';'])
    eval(['M1(2,:)=[' names{m,1} '_' Annotators{2} '];'])
%     M1(isnan(M1)) = 43 ; 
    FK1 = confusionmat(M1(1,:),M1(2,:))  ;

    [k,sek]=kappa(FK1,2); % k=kappa, sek= kappa error
    
    eval(['k_' num2str(pat) '(n,1)=k;'])
    
end
    
    
      clearvars -except p Neonate Annotators savefolder k_*
  
end