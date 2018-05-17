function     RemoveNaNoverallMatrices(Neonate,win,saving,Loadfolder,savefolder)

disp('Removing NAN in matrices')

    for j=1:length(win)
        Neonate;
        disp(['working on window: ' num2str(win(1,j))]);

        folderF=dir([Loadfolder 'Feature_Matrix_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']);
        folderA=dir([Loadfolder 'Annotation_Matrix_*' '_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']);
        folderW=dir([Loadfolder 'sampleWeight_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']);
               
        
        
        for k=1:length(folderF)
            FMatrix=cell2mat(struct2cell((load([Loadfolder folderF(k,1).name],'Matrix')))); 
            disp([folderF(k,1).name])
            [row, col]=find(all(isnan(FMatrix),1));
            FMatrix(:,col)=[];
            Matrix=FMatrix;   
            save ([savefolder folderF(k,1).name],'Matrix')
            Matrix=[];
                for i=1:length(folderA)
                    disp([folderA(i,1).name])
                    AMatrix=cell2mat(struct2cell((load([Loadfolder folderA(i,1).name],'Matrix')))); 
                    AMatrix(:,col)=[];
                    Matrix=AMatrix;   
                    save ([savefolder folderA(i,1).name],'Matrix')            
                    Matrix=[];            
                end
                for p=1:length(folderW)        
                    disp([folderW(p,1).name])
                    WMatrix=cell2mat(struct2cell((load([Loadfolder folderW(p,1).name],'Matrix')))); 
                    WMatrix(:,col)=[];
                    Matrix=WMatrix;   
                    save ([savefolder folderW(p,1).name],'Matrix') 
                    Matrix=[];
                end
        end        
    end

end
