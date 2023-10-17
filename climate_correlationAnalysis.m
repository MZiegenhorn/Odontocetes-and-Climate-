%created by MAZ on 10/11/2022 

clear all; close all; clc

%grab climate indices values and times
load('allClimate.mat')

%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%
%load in variable of interest
infile = '/Volumes/MAZ_keep/ch3/correlation';
outfile = '/Users/maziegenhorn/Documents/thesis_chapters/ch3/correlation';
varname = 'chla';
site = 'all'; %hawaii, manawai, or all
lag = 0;
pval = 0.05; %what do you want to set significance level less than for maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isdir(outfile)
    mkdir(outfile)
end

load([infile,'/',varname,'month.mat']);

dataexp = cat(3,monthvar{:});

%setup climate variables for correlation purposes
ensoval = [];
pdoval = [];
npgoval = [];

%for chla,have 2007 data which we don't for others, so remove to start with
%july 2008 (as for hycom)
jul2008 = datenum(2008,7,0,0,0,0);
rmvmonth = find(unmonth<jul2008);
unmonth(rmvmonth) = [];
dataexp(:,:,rmvmonth) = [];

%deal with site geolimits
if strcmp(site,'all')
    latlims = [19 29]; lonlims =  [-177 -155];
    siteloc = [27.7326 -175.6036];
elseif strcmp(site, 'Hawaii')
    latlims = [19 23]; lonlims = [-160 -155.5];
    siteloc = [0,0];
elseif strcmp(site, 'manawai')
    latlims = [24 30]; lonlims = [-180 -170];
    siteloc = [27.7326 -175.6036];
end


for iu = 1:size(unmonth,1)
    dvmonth = datevec(unmonth(iu));
    year = dvmonth(1); eyear = find(enso(:,1) == year);
    month = dvmonth(2);
    ensoval(iu) = enso(eyear,month+1);
    
    %pdo
    pyear = find(pdo(:,1) == year);
    pdoval(iu) = pdo(pyear,month+1);
    
    %npgo
    nrow = find(npgo(:,1) == year & npgo(:,2) == month);
    npgoval(iu) = npgo(nrow,3);
end

ensoLoc = [];
pdoLoc = [];
npgoLoc = [];
ensoP = [];
pdoP = [];
npgoP = [];

if ~strcmp(varname,'chla')
    
    for ia = 1:size(dataexp,1)
        %for each lat lon pair, look at correlation
        lonuse = lons(ia);
        
        for ib = 1:size(dataexp,2)
            latuse = lats(ib);
            %data is stored in lon, lat (i.e. there are 276 longitude points)
            datause = dataexp(ia,ib,:);
            data = reshape(datause,1,[]);
            
            %test correlation between this value and value of climate variabls in
            %corresponding month
            [Re,Pe] = corrcoef(data,ensoval);
            [Rp, Pp] = corrcoef(data,pdoval);
            [Rn, Pn] = corrcoef(data,npgoval);
            
            ensoLoc(ia,ib) = Re(2);
            pdoLoc(ia,ib) = Rp(2);
            npgoLoc(ia,ib) = Rn(2);
            ensoP (ia,ib)= Pe(2);
            pdoP(ia,ib) = Pp(2);
            npgoP(ia,ib) = Pn(2);
        end
    end
    
    %prune out anything with p>0.05
    ensoLoc(ensoP>pval) = NaN;
    pdoLoc(pdoP>pval) = NaN;
    npgoLoc(npgoP>pval) = NaN;
    
    
    %plot correlation values on our map grid
    ln = repmat(lons',1,length(lats));
    lt = repmat(lats,1,length(lons))';
    
    figure
    
    subplot(2,2,1)
    pcolor(ln,lt,ensoLoc)
    hold on
    plot(siteloc(2)+360,siteloc(1),'.k','MarkerSize',40)
    shading interp
    %set the colorbar axis to just include values we want
    cmocean('curl')
    caxis([-0.6 0.6])
    colorbar
    axis tight
    title('ENSO Corr')
    
    subplot(2,2,2)
    pcolor(ln,lt,pdoLoc)
    hold on
    plot(siteloc(2)+360,siteloc(1),'.k','MarkerSize',40)
    shading interp
    %set the colorbar axis to just include values we want
    cmocean('curl')
    caxis([-0.6 0.6])
    colorbar
    axis tight
    title('PDO Corr')
    
    subplot(2,2,3)
    pcolor(ln,lt,npgoLoc)
    hold on
    plot(siteloc(2)+360,siteloc(1),'.k','MarkerSize',40)
    shading interp
    %set the colorbar axis to just include values we want
    cmocean('curl')
    caxis([-0.6 0.6])
    colorbar
    axis tight
    title('NPGO Corr')
    
else
    
    if lag > 0
        %use the lag to shift our data
        %if there's a 6 month lag, we want to know what chla will be 6
        %months from now, so remove first lag number of data points and go
        %from there
        dataexp(:,:,1:lag) = [];
        %remove last 6 values from enso, pdo, npgo values
        ensoval = ensoval(1:end-lag);
        pdoval = pdoval(1:end-lag);
        npgoval = npgoval(1:end-lag);
    end
    
    for ia = 1:size(dataexp,1)
        %for each lat lon pair, look at correlation
        lonuse = lons(ia);
        latuse = lats(ia);
        %data is stored in lon, lat (i.e. there are 276 longitude points)
        datause = dataexp(ia,1,:);
        data = reshape(datause,1,[]);
        
        %test correlation between this value and value of climate variabls in
        %corresponding month
        [Re, Pe] = corrcoef(data,ensoval,'Rows','complete');
        [Rp, Pp] = corrcoef(data,pdoval,'Rows','complete');
        [Rn, Pn] = corrcoef(data,npgoval,'Rows','complete');
        
        ensoLoc(ia) = Re(2);
        pdoLoc(ia) = Rp(2);
        npgoLoc(ia) = Rn(2);
        ensoP(ia) = Pe(2);
        pdoP(ia) = Pp(2);
        npgoP(ia) = Pn(2);
    end
    
    %get rid of non-significant correlations
    ensoLoc(ensoP>pval) = NaN;
    pdoLoc(pdoP>pval) = NaN;
    npgoLoc(npgoP>pval) = NaN;
    
    figure
    
    subplot(2,2,1)
    geoscatter(lats,lons,5,ensoLoc,'filled')
    hold on
    geoscatter(siteloc(1),siteloc(2),80,1,'filled')
    geolimits(latlims, lonlims)
    %set the colorbar axis to just include values we want
    cmocean('curl')
    caxis([-0.6 0.6])
    colorbar
    title('ENSO Corr')
    
    %pdo
    subplot(2,2,2);
    geoscatter(lats,lons,10,pdoLoc,'filled')
    hold on
    geoscatter(siteloc(1),siteloc(2),80,1,'filled')
    geolimits(latlims, lonlims)
    %set the colorbar axis to just include values we want
    cmocean('curl')
    caxis([-0.6 0.6])
    colorbar
    title('PDO Corr')
    
    %npgo
    subplot(2,2,3);
    geoscatter(lats,lons,10,npgoLoc,'filled')
    hold on
    geoscatter(siteloc(1),siteloc(2),80,1,'filled')
    geolimits(latlims, lonlims)
    %set the colorbar axis to just include values we want
    cmocean('curl')
    caxis([-0.6 0.6])
    colorbar
    title('NPGO Corr')
    
end

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.8, 0.9]);
sgtitle(varname,'Interpreter', 'none')
if lag>0
    print([outfile,'/',varname,'ClimateCorr_lag',num2str(lag),'_',site],'-dpng')
    
    save([outfile,'/',varname,'_lag',num2str(lag),'_',site,'_correlations.mat'],'ensoLoc','pdoLoc','npgoLoc',...
        'ensoval','pdoval','npgoval','lons','lats')
else
    print([outfile,'/',varname,'ClimateCorr_',site],'-dpng')
    
    save([outfile,'/',varname,'_',site,'_correlations.mat'],'ensoLoc','pdoLoc','npgoLoc',...
        'ensoval','pdoval','npgoval','lons','lats')
end



