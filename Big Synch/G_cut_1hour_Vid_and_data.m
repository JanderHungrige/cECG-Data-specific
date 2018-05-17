clc
clear

ffmpegfolder=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Quick Annotator');

pat=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
pat=[1,2,3,5,7,8,9,10,11,12,13,14];
pat=2;


for i=1:length(pat)
  disp(['splitting data for patient ' num2str(pat(1,i))])
        Neonate=pat(i);
        cd(ffmpegfolder)
    split_data_in_1hour(Neonate)
    
end % for each patient
