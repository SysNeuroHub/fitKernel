[saveServer, rootFolder] = getReady();

loadNames_all = [];monthdates_all = [];
for aa = 1:2
    switch aa
        case 1
            animal =  'hugo';%'ollie';% % %'andy';%
        case 2
            animal = 'ollie';
    end
    for yyy = 1:3
        switch yyy
            case 1
                year = '2021'; %hugo
            case 2
                year = '2022'; %hugo
            case 3
                year = '2023'; %hugo, ollie
        end

        [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

        monthdates = [];
        for ii = 1:numel(dates)
            monthdates{ii} =  [months{ii} dates{ii}];
        end

        monthdates_all = [monthdates_all unique(monthdates)];
        loadNames_all = [loadNames_all loadNames];
    end
end