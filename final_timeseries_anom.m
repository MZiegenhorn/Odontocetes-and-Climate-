%created by MAZ on 6/9/22 to explore data timeseries for final model
%variables only
%https://www.mathworks.com/matlabcentral/fileexchange/68969-colorpalette


clear all
close all
clc

%%%%%%%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%
site = 'Hawaii';
infile = '/Users/maziegenhorn/Documents/SIO/thesis_chapters/ch3/manuscript/final/code_and_data/Manawai_allclasses_counts_finalData.mat';
outfolder = ['/Users/maziegenhorn/Documents/thesis_chapters/ch3/finalTS/',site];
class = {'1','2','3','4','5_6','7_8','9','10'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isdir(outfolder)
    mkdir(outfolder)
end

ksvalhisto = [];
ksvaldata = [];

%open matfile
load(infile)

varnames = daydata.Properties.VariableNames;

%if the column for daycounts is not at start, move it
daycountsloc = find(strcmp(varnames,'daycounts'));
if daycountsloc > 2
    daydata = [daydata(:,daycountsloc),daydata(:,1:daycountsloc-1),daydata(:,daycountsloc+1:end)];
end

dayarray = table2array(daydata);
varnames = daydata.Properties.VariableNames;

timeuse = daydata;
timearray = table2array(timeuse);


%run through for each type
for ic = 1:size(class,2)
    %set up variable names
    varsuse = {'pdo','enso','npgo','ssht_dep0','salt_dep0','tempt_dep0'};
    %add variable for current class
    varsuse = [{'day'},{['class',class{ic}]},varsuse];

    %truncate daydata
    varind = [];
    for iv = 1:size(varsuse,2)
        varind(iv) = find(strcmp(varnames,varsuse{iv}));
    end

    dayuse = timeuse(:,varind);
    dayarrayfinal = table2array(dayuse);
    varsforplot = dayuse.Properties.VariableNames;

    if strcmp(site,'manawai')
        %split at start of man08
        dayusek = dayuse;
        daystr = datevec(dayusek.day);
        man08start = min(find(daystr(:,1) == 2014));
        stends = [1,man08start-1;man08start,length(daystr)];
        for ik = 1:2
            dayuse = dayusek(stends(ik,1):stends(ik,2),:);
            dayarrayfinal = table2array(dayuse);
            %get a full dayuse for plotting- we can use nans where necessary
            alldays = min(dayuse.day):datenum(0,0,1,0,0,0):max(dayuse.day);
            [~,ia] = intersect(alldays,dayuse.day);
            %get boxes of no effort setup
            ids = ones(size(alldays));
            ids(ia) = 0;
            offdays = alldays(logical(ids));
            if ~isempty(offdays)
                daydiff = [offdays(2:end),offdays(end)] - offdays;
                %find large gaps
                gaps = find(daydiff > 10); %find any gaps longer than a week
                gaps = [gaps,size(offdays,2)]; %make sure to get end in there
            end
            %actual plotting
            figure
            step = 1;
            nrow = size(dayuse,2)-1;
            for ic2 = 2:size(dayuse,2)
                subplot(nrow,1,step)

                ylab = varsforplot{ic2};

                %get a full dayuse for plotting- we can use nans where necessary
                alldays = min(dayuse.day):datenum(0,0,1,0,0,0):max(dayuse.day);
                fullvar = nan(size(alldays));
                %get anomaly for plotting
                varuse = dayarrayfinal(:,ic2);
                %setup array for getting means across both
                arrayformean = table2array(dayusek);
                varuse1 = arrayformean(:,ic2);
                %find the spots where variable should be populated
                [~,ia] = intersect(alldays,dayuse.day);

                %if we're looking at climatic variable, put a line at 0
                if contains(varsforplot{ic2},'enso')||contains(varsforplot{ic2},'pdo')...
                        ||contains(varsforplot{ic2},'npgo')
                    meanval = 0;
                elseif contains(varsforplot{ic2},'class')
                    %if we're looking at class, get weekly values and do
                    %anomalies
                    %get means to be same across both plots
                    allbins = min(dayusek.day):datenum(0,0,7,0,0,0):max(dayusek.day)+datenum(0,0,7,0,0,0);
                    [~,~,iduse1] = histcounts(dayusek.day,allbins);
                    varuse1 = accumarray(iduse1,arrayformean(:,ic2));
                    meanval = mean(varuse1(unique(iduse1)));
                    %get actual bins for plots
                    weekbins = min(dayuse.day):datenum(0,0,7,0,0,0):max(dayuse.day)+datenum(0,0,7,0,0,0);
                    [binnedays,~,iduse] = histcounts(dayuse.day,weekbins);
                    varuse = accumarray(iduse,dayarrayfinal(:,ic2));
                    %                         meanval = mean(varuse(unique(iduse1)));
                    alldays = weekbins(1:end-1);
                    fullvar = nan(size(alldays));
                    ia = 1:size(alldays,2);
                else
                    %otherwise, put the line at the mean
                    meanval = mean(varuse1);
                end
                negvar = fullvar;
                negvar(ia(varuse<meanval)) = varuse(varuse<meanval);
                posvar = fullvar;
                posvar(ia(varuse>meanval)) = varuse(varuse>meanval);

                %handle color
                if contains(varsforplot{ic2},'class')
                    col1 = [92/255,188/255,99/255];%'fern'
                    col2 = [150/255,250/255,169/255];%'lighter... fern...'
                elseif contains(varsforplot{ic2},'enso')||contains(varsforplot{ic2},'pdo')...
                        ||contains(varsforplot{ic2},'npgo')
                    col1 = [98/255,197/255,218/255]; %sky
                    col2 = [40/255,30/255,93/255];%indigo
                elseif contains(varsforplot{ic2},'sal')
                    col1 = [240/255,128/255,128/255];%'orange red'
                    col2 = [255/255,218/255,185/255];%'orange-pink'
                elseif contains(varsforplot{ic2},'temp')
                    col1 = [236/255,151/255,6/255];%'honey'
                    col2 = [253/255,238/255,135/255];%'daffodil'
                elseif contains(varsforplot{ic2},'ssh')
                    col1 = [204/255,0/255,0/255];%'darkred'
                    col2 = [255/255,153/255,102/255];%'red'
                end

                %plot the variable against time-separate pos and negative
                if contains(varsforplot{ic2},'class')
                    area(alldays,varuse,0,'FaceColor',col1)
                    hold on
                else
                    area(alldays,posvar,meanval,'FaceColor',col1)
                    hold on
                    area(alldays,negvar,meanval,'FaceColor',col2)
                end
                %add patches for empty data
                if ~isempty(offdays)
                    step2 = 1;
                    for ig = 1:size(gaps,2)
                        patch([offdays(step2) offdays(step2) offdays(gaps(ig)) offdays(gaps(ig))],...
                            [min(varuse1) max(varuse1) max(varuse1) min(varuse1)],[0.92 0.92 0.92],'EdgeColor',[0.9 0.9 0.9])
                        step2 = gaps(ig) + 1;
                    end
                end

                hold off

                ylabel(ylab)
                xlim([min(dayuse.day) max(dayuse.day)])

                if max(varuse1)~=min(varuse1)
                    ylim([min(varuse1),max(varuse1)])
                end
                %set x-tick locations to be at start of each year
                lowyr = datevec(min(dayuse.day));
                highyr = datevec(max(dayuse.day));
                years = datenum(lowyr(1),0,0,0,0,0):datenum(1,0,0,0,0,0):datenum(highyr(1)+1,0,0,0,0,0);
                %                     years = [years,years+datenum(0,6,0,0,0,0)];
                %                     xticks(sort(years))
                %                     ax = gca; ax.LineWidth =0.8; % Make xtick thicker
                if ic2 ~= size(dayuse,2)
                    set(gca,'XTickLabels',[])
                end
                step = step + 1;

            end
            %save our current stuff, reset step
            datetick('x','mm/dd/yy','keepticks')
            xlim([min(dayuse.day) max(dayuse.day)])
            xlabel('Time')
            sgtitle(['Timeseries for class',class{ic}])

            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.01, 0.7, 1]);
            %           save image
            plotname = [outfolder,'\',site,'_',num2str(ik),'class',class{ic},'allTS_final'];
            print(plotname,'-dpng')
        end
    else

        %get a full dayuse for plotting- we can use nans where necessary
        alldays = min(timeuse.day):datenum(0,0,1,0,0,0):max(timeuse.day);
        [~,ia] = intersect(alldays,dayuse.day);
        %            get boxes of no effort setup
        ids = ones(size(alldays));
        ids(ia) = 0;
        offdays = alldays(logical(ids));
        daydiff = [offdays(2:end),offdays(end)] - offdays;
        %find large gaps
        gaps = find(daydiff > 10); %find any gaps longer than a week
        gaps = [gaps,size(offdays,2)]; %make sure to get end in there

        %actual plotting
        figure
        step = 1;
        nrow = size(dayuse,2)-1;
        for ic2 = 2:size(dayuse,2)
            subplot(nrow,1,step)

            ylab = varsforplot{ic2};

            %get a full dayuse for plotting- we can use nans where necessary
            alldays = min(timeuse.day):datenum(0,0,1,0,0,0):max(timeuse.day);
            fullvar = nan(size(alldays));
            %get anomaly for plotting
            varuse = dayarrayfinal(:,ic2);
            %find the spots where variable should be populated
            [~,ia] = intersect(alldays,dayuse.day);

            %if we're looking at climatic variable, put a line at 0
            if contains(varsforplot{ic2},'enso')||contains(varsforplot{ic2},'pdo')...
                    ||contains(varsforplot{ic2},'npgo')
                meanval = 0;
            elseif contains(varsforplot{ic2},'class')
                %if we're looking at class, get weekly values and do
                %anomalies
                weekbins = min(dayuse.day):datenum(0,0,7,0,0,0):max(dayuse.day)+datenum(0,0,7,0,0,0);
                [binnedays,~,iduse] = histcounts(dayuse.day,weekbins);
                varuse = accumarray(iduse,dayarrayfinal(:,ic2));
                meanval = mean(varuse(unique(iduse)));
                alldays = weekbins(1:end-1);
                fullvar = nan(size(alldays));
                ia = 1:size(alldays,2);
            else
                %otherwise, put the line at the mean
                meanval = mean(varuse);
            end
            negvar = fullvar;
            negvar(ia(varuse<meanval)) = varuse(varuse<meanval);
            posvar = fullvar;
            posvar(ia(varuse>meanval)) = varuse(varuse>meanval);

            %handle color
            if contains(varsforplot{ic2},'class')
                col1 = [92/255,188/255,99/255];%'fern'
                col2 = [150/255,250/255,169/255];%'lighter... fern...'
            elseif contains(varsforplot{ic2},'enso')||contains(varsforplot{ic2},'pdo')...
                    ||contains(varsforplot{ic2},'npgo')
                col1 = [98/255,197/255,218/255]; %sky
                col2 = [40/255,30/255,93/255];%indigo
            elseif contains(varsforplot{ic2},'sal')
                col1 = [240/255,128/255,128/255];%'orange red'
                col2 = [255/255,218/255,185/255];%'orange-pink'
            elseif contains(varsforplot{ic2},'temp')
                col1 = [236/255,151/255,6/255];%'honey'
                col2 = [253/255,238/255,135/255];%'daffodil'
            elseif contains(varsforplot{ic2},'ssh')
                col1 = [204/255,0/255,0/255];%'darkred'
                col2 = [255/255,153/255,102/255];%'red'
            end
            %plot the variable against time-separate pos and negative
            if contains(varsforplot{ic2},'class')
                area(alldays,varuse,0,'FaceColor',col1)
                hold on
            else
                area(alldays,posvar,meanval,'FaceColor',col1)
                hold on
                area(alldays,negvar,meanval,'FaceColor',col2)
            end
            %add patches for empty data
            step2 = 1;
            for ig = 1:size(gaps,2)
                patch([offdays(step2) offdays(step2) offdays(gaps(ig)) offdays(gaps(ig))],...
                    [min(varuse) max(varuse) max(varuse) min(varuse)],[0.92 0.92 0.92],'EdgeColor',[0.9 0.9 0.9])
                step2 = gaps(ig) + 1;
            end

            hold off
            %                 plot(alldays,zeros(size(alldays)),'b')

            ylabel(ylab)
            xlim([min(dayuse.day) max(dayuse.day)])
            if max(varuse)~=min(varuse)
                ylim([min(varuse),max(varuse)])
            end
            %set x-tick locations to be at start of each year
            lowyr = datevec(min(dayuse.day));
            highyr = datevec(max(dayuse.day));
            years = datenum(lowyr(1),0,0,0,0,0):datenum(1,0,0,0,0,0):datenum(highyr(1)+1,0,0,0,0,0);
            years = [years,years+datenum(0,6,0,0,0,0)];
            xticks(sort(years))
            if ic2 ~= size(dayuse,2)
                set(gca,'XTickLabels',[])
            end
            step = step + 1;

        end
        %save our current stuff, reset step
        datetick('x','mm/dd/yy','keepticks')
        xlim([min(dayuse.day) max(dayuse.day)])
        xlabel('Time')
        sgtitle(['Timeseries for class',class{ic}])

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.01, 0.7, 1]);
        %           save image
        plotname = [outfolder,'\',site,'_class',class{ic},'allTS_final'];
        print(plotname,'-dpng')
        %     end
    end
end

