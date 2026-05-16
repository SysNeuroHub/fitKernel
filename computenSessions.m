nUnits = [];
nSession = [];
for aa = 1:2
    switch aa
        case 1
            animal = 'hugo'; yidx=1:3;
        case 2
            animal =  'ollie'; yidx = 3;
    end
    for yy = yidx
        switch yy
            case 1
                year = '2021';
            case 2
                year = '2022';
            case 3
                year = '2023';
        end
        [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);
        nUnits(aa,yy) = numel(loadNames);
        

        thisDate = [];
        for idata = 1:nUnits(aa,yy)
            thisDate{idata} = [months{idata} '_' dates{idata}];
        end
        nSession(aa,yy) = numel(unique(thisDate));
    end
end