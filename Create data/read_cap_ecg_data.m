%Aline's Code, reading reference respiration and capacitive data.
function [referenceRespi,SEN]= read_cap_ecg_data(directory)
%directory='C:\Projects\neonatal\data\March2012\neonatal data';
cd(directory)

fileNumber=0;

% Get size of the file
        fid = fopen(['Ref ECG1_0000000',num2str(fileNumber),'.txt'], 'r');
        temp = fread(fid,inf,'*char')';
        a=8; sizeFile=NaN;
        while isnan(sizeFile)
            sizeFile=str2double(temp(end-a:end));
            a=a-1;
        end
        fclose(fid);

sampleCounter=0;
numberNewSamples=sizeFile;



% reference ECG
            mref = eval(['memmapfile(''Ref ECG1_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
            xref = mref.data(3*sampleCounter+1:3*sampleCounter+3*numberNewSamples); %binary data
            LengthData = length(xref); % in bytes
            %reading the binary files then puts in ECG, ref data
            ECG1 = double(typecast(xref(1:3:LengthData),'int8')')*2^16+double(xref(2:3:LengthData)')*256+double(xref(3:3:LengthData)');
            % reference respi
            mref = eval(['memmapfile(''Ref Resp_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
            xref = mref.data(3*sampleCounter+1:3*sampleCounter+3*numberNewSamples);
            LengthData = length(xref); % in bytes
            referenceRespi = double(typecast(xref(1:3:LengthData),'int8')')*2^16+double(xref(2:3:LengthData)')*256+double(xref(3:3:LengthData)');
            % capacitive sensors
            for sens=1:8
                m = eval(['memmapfile(''Sensor ',num2str(sens),'_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
                x = m.data(3*sampleCounter+1:3*sampleCounter+3*numberNewSamples);
                eval(['SEN',num2str(sens),' = double(typecast(x(1:3:LengthData),''int8'')'')*2^16+double(x(2:3:LengthData)'')*256+double(x(3:3:LengthData)'');']);
            end




Siglist={'SEN1','SEN2','SEN3','SEN4','SEN5','SEN6','SEN7','SEN8'};
for i=1:8
    SEN{i}=eval(Siglist{i});
end
