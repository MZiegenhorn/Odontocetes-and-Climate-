%created by MAZ on 5/3/22 to add climate indices to habitat modelling spreadsheet
 
clear all
close all
clc

%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

site = 'manawai'; %site for focus
indata = ['E:\ch3\modelling\',site]; %folder of input data

%hycom
hycom = 'E:\ch3\HYCOM'; %where hycom data is stored
%depths to pull HYCOM data from
usedepthsnew = [-0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global depths


allfiles = dir(fullfile(indata,['*binsperhr','*habMod.mat']));

for ia = 1:size(allfiles,1)
    %open file
    load(fullfile(allfiles(ia).folder,allfiles(ia).name))
    usedepths = usedepthsnew;
    days = daydata.day;
    
    %%hycom-just one value per day, so not local time
    hyfile = fullfile(hycom,[site,'_allVars']);
    load(hyfile)
    [~,depthinds] = intersect(depths,usedepths);
    ssht = struct2table(get_hycomvar(ssh,'ssh',1,days,time));
    salt = struct2table(get_hycomvar(sal,'sal',depthinds,days,time));
    tempt = struct2table(get_hycomvar(temp,'temp',depthinds,days,time));
    
    %modify these in our structure
    %remove current hycom fields
    varnames = daydata.Properties.VariableNames;
    rmvfields = contains(varnames,'_dep');
    daydata(:,rmvfields) = [];
    %replace with new hycom variables
    hytable = table(ssht,salt,tempt);
    %use splitvars to get variables out
    hytable = splitvars(hytable);
    daydata = [daydata,hytable];
    
    %save the data
    filename = strrep(allfiles(ia).name,'.mat','.csv');
    saven = fullfile(indata,filename);
    writetable(daydata,saven,'Delimiter',',')
    %save as mat file for other plotting
    saven2 = strrep(saven,'.csv','.mat');
    save(saven2,'daydata','usedepths','sitelat','sitelon','nanuse','nonnan')
    
    disp(['Done with ',allfiles(ia).name])
end


function outvar = get_hycomvar(var,varname,depthinds,day,time)
%find right data, and assign it correctly
global depths
outvar = struct();
depthinds = flipud(depthinds);
nanflag = [];
for idh = 1:length(day)
    %find the day in our data
    usedayh = find(time==day(idh));
    depthsv = [];
    %extract our variable at that point-run through each to make variables
    for idepth = 1:size(depthinds,1)
        depthn = ['dep',num2str(abs(depths(depthinds(idepth))))];
        if isempty(usedayh)
            disp(['WARNING: no ',varname,' data found for ',datestr(day(idh))])
            outvar.(depthn)(idh,:) = NaN;
            nanflag(idh) = 1;
        else
            outvar.(depthn)(idh,:) = var(usedayh,depthinds(idepth));
            nanflag(idh) = 0;
        end
        depthsv{idepth} = depthn;
    end
end

%interpolate for temporal nans
nanuse = find(nanflag == 1);
nonnan = find(nanflag == 0);
if ~isempty(nanuse)
    disp('Averaging for temporal NaNs')
    for idep = 1:size(depthsv,2)
        outtemp = interp1(day(nonnan),outvar.(depthsv{idep})(nonnan),day,'linear');
        if size(outtemp,1) == 1
            outvar.(depthsv{idep}) = outtemp';
        else
            outvar.(depthsv{idep}) = outtemp;
        end
    end
end
%
disp(['Done with variable ',varname])
end

