function [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder)
%[loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder)

foldersByYear = dir(rootFolder);
foldersByYear = foldersByYear(3:end);
foldersByYear = foldersByYear([foldersByYear.isdir]);
thisFolder = strcmp({foldersByYear.name}, year);
if sum(thisFolder)==0
    error('no data found in this year');
end
foldersByYear = foldersByYear(thisFolder);

thisFolder = fullfile(foldersByYear.folder,foldersByYear.name, 'cuesaccade_data');
foldersByMonth = dir(thisFolder);
foldersByMonth = foldersByMonth(3:end);
foldersByMonth = foldersByMonth([foldersByMonth.isdir]);

monthNames = {foldersByMonth.name};
dateNames = cell(1,length(monthNames));
for im = 1:length(monthNames)
    foldersByDate = dir(fullfile(thisFolder, monthNames{im}));
    foldersByDate = foldersByDate(3:end);
    foldersByDate = foldersByDate([foldersByDate.isdir]);
    dateNames{im} = {foldersByDate.name};
end

idata = 1;
loadNames = [];
months = [];
dates = [];
channels = [];
for im = 1:length(monthNames)
    for id = 1:length(dateNames{im})
        date = [monthNames{im} '/' dateNames{im}{id}];
        
        loadFolder = fullfile(thisFolder, date, 'saved_oephysdata');
        
        fileNames = {dir(fullfile(loadFolder, [animal '_oephysdata_*.mat'])).name};
        
        for ich = 1:length(fileNames)
            loadName = fullfile(loadFolder, fileNames{ich});
            
            if exist(loadName, 'file')
                loadNames{idata} = loadName;
                months{idata} = monthNames{im};
                dates{idata} = dateNames{im}{id};
                chName_c =  textscan(fileNames{ich}, [animal '_oephysdata_ch%d.mat']);
                channels{idata} = chName_c{1};
                idata = idata + 1;
            end
        end
    end
end