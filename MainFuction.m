% p = parpool('local',1);
tic;
SubFolder = GetFolders('.\TargetSpectrum');
eqname={'EW.txt'};%name of earthquake file
for ii=1:length(SubFolder)
    disp(ii);
    temp_string=strsplit(SubFolder{1,ii},'\');
    cur_foldername=temp_string{length(temp_string)};
    this_station=char(cur_foldername);
    dsname   = 'PredictedGM.txt';                 % name of spectrum file (2 columns
    eqfolder = '.\groundmotion';   % directory with the accelerograms 输出地震动目录
    ArtifQuakeLetII(SubFolder{1,ii},dsname,eqfolder,eqname{1,1},this_station);
end
toc;
% delete(p);



