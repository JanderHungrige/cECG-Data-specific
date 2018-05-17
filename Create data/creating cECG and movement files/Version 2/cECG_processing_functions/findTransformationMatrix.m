function [vcg2ecgMatrix,transformationMatrix]=findTransformationMatrix(typeOfArray, channelsToUse)

%spatial combination of signals based on their electrode positions. This
%routine gives the same matrix everytime, so in the main file this fixed
%matrix is used. This script nevertheless contains the derivation of the
%fixed matrix.

switch typeOfArray
    case 'heptagone'
        % Perfect neonatal heptagone
        electrodes(1,:)=[15,31.15];
        electrodes(2,:)=[33.7,7.69];
        electrodes(3,:)=[27.03,-21.55];
        electrodes(4,:)=[0,-34.57];
        electrodes(5,:)=[-27.03,-21.55];
        electrodes(6,:)=[-33.7,7.69];
        electrodes(7,:)=[-15,31.15];
        electrodes(8,:)=[0,0];
        
    case 'neonatalV1'
        % Real neonatal configuration
        electrodes(1,:)=[34,28.5];
        electrodes(2,:)=[42.5,-3];
        electrodes(3,:)=[28.5,-30.5];
        electrodes(4,:)=[-3,-41];
        electrodes(5,:)=[-32,-23];
        electrodes(6,:)=[-30,15.5];
        electrodes(7,:)=[-1,34];
        electrodes(8,:)=[0,0];
        
    case 'neonatalV2'
        % Perfect neonatal heptagone shifted 180 degrees
%         electrodes(5,:)=[15,31.15];
%         electrodes(6,:)=[33.7,7.69];
%         electrodes(7,:)=[27.03,-21.55];
%         electrodes(1,:)=[0,-34.57];
%         electrodes(2,:)=[-27.03,-21.55];
%         electrodes(3,:)=[-33.7,7.69];
%         electrodes(4,:)=[-15,31.15];
%         electrodes(8,:)=[0,0];
        
        electrodes(4,:)=[15,31.15];
        electrodes(5,:)=[33.7,7.69];
        electrodes(6,:)=[27.03,-21.55];
        electrodes(7,:)=[0,-34.57];
        electrodes(1,:)=[-27.03,-21.55];
        electrodes(2,:)=[-33.7,7.69];
        electrodes(3,:)=[-15,31.15];
        electrodes(8,:)=[0,0];
    case 'adult'
        % Adult config
        electrodes(1,:)=[39,72];
        electrodes(2,:)=[81,15];
        electrodes(3,:)=[60,-54];
        electrodes(4,:)=[1,-78];
        electrodes(5,:)=[-67,-45];
        electrodes(6,:)=[-81,26];
        electrodes(7,:)=[-35,80];
        electrodes(8,:)=[0,0];   
end

vcg2ecgMatrix=electrodes(channelsToUse,:)-repmat(mean(electrodes(channelsToUse,:),1),[numel(channelsToUse) 1]); 

% angles = atan(vcg2ecgMatrix(:,2)./vcg2ecgMatrix(:,1));
% angles(6:8) = angles(6:8)+pi;
% angles(4:5) = angles(4:5)-pi;
% transfoMat = [cos(angles)' ; sin(angles)'];

vcg2ecgMatrix=vcg2ecgMatrix./repmat(sum(abs(vcg2ecgMatrix),1),[size(vcg2ecgMatrix,1) 1]); %normalization %Why is this needed? Why not doing it with the transformationMatrix? 
transformationMatrix=inv(vcg2ecgMatrix'*vcg2ecgMatrix)*vcg2ecgMatrix';
%transformationMatrix =  vcg2ecgMatrix'; %Why not using this instead of
%previous line?
