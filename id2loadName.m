function loadName = id2loadName(thisid, rootFolder)
% thisid = 
% 'hugo/2021/03March/16/13'
%
% datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
% thisid = [animal '/' year '/' datech];


% loadName = 
% '/mnt/MBI/Monash
% Data/Joanita/2021/cuesaccade_data/03March/16/saved_oephysdata/hugo_oephysdata_ch13.mat"


parts = strsplit(thisid, filesep);
animal = parts{1};
year = parts{2};
month = parts{3};
date = parts{4};
ch = parts{5};

loadFolder = fullfile(rootFolder, year, 'cuesaccade_data',month,date, 'saved_oephysdata');
fileName = [animal '_oephysdata_ch' ch '.mat'];
loadName = fullfile(loadFolder, fileName);
