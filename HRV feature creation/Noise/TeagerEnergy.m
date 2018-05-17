%teager energy

Neonate=3;
cd('C:\Users\310122653\Documents\PhD\EEG_MMC Data\Matlab')
[sheet,lead,FS]=callingSheetsLeadsandFS(Neonate); % calling which sheet has to be chosen for the annotation, which lead for the ECG and which FS results from that

load(['C:\Users\310122653\Documents\PhD\InnerSense Data\Data files\All files\Baby' num2str(Neonate) '\sfull']); %Full set

%%% load data
ECG_raw=s.Data{1,lead};
QRS=ECG_raw(1,150:300);
noise=ECG_raw(1,3.714e6:3.72E6);
plotting=1;
%teagerenergy
cd('C:\Users\310122653\Documents\PhD\EEG_MMC Data\Matlab\HRV anlaysis\Noise')

[ey,ex]=energyop(noise,plotting);
