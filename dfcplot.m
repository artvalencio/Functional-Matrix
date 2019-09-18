function [links,lags,ds1]=dfcplot(eegdata,chanlocs,shift,windowlen,step,dfctype,threshold)
%DFCPLOT Makes the grid with DFC at the times selected
%Requires EEGLAB and Praneeth Namburi's Connected Topoplot
%-----------------------------------------
%Inputs:
%- eegdata:         the EEG data, given in EEGLAB format
%- chanlocs:        the electrode locations, in the format provided by EEGLAB
%                   (after importing data check EEG.chanlocs variable)
%- shift:           the maximum shift allowed for correlation analysis considering 
%                   time-delay between locations. Given in time-series points
%                   units. Select shift=0 for disregard time-delay effects
%- windowlen:       the interval (window length) on which the correlation or information is 
%                   calculated
%- step:            step (jump) given between calculation of correlation in one window
%                   and on the other. Select step=1 for better precision,
%                   step>1 for faster computation
%- dfctype:         'correl' for calculating with Pearson correlation, 'spearman' for Spearman
%                   correlation, 'kendall' for Kendall correlation, 'info' for
%                   information-theoretical measures (CaMI, Mutual Information,Transfer
%                   Entropy,and Directionality)
%- threshold:       defines the threshold beyond which a correlation/information
%                   measure is considered a valid link.
%                   e.g. threshold=0.75 for dfctypes 'correl','spearman' or
%                   'kendall' means that correlations above 0.75 or below -0.75
%                   are interpreted as establishing a link
%                   for dfctype='info', threshold must be given as a 4-element vector:
%                   [thr_cami thr_mi thr_te thr_dir], with each element corresponding 
%                   to the threshold value of a respective information measure
%-----------------------------------------------------------------------------------
%Outputs
%Figures with the topoplots showing the connectivities in each selected
%time interval
%-----------------------------------------------------------------------------------
%Usage example
%[links,lags,ds1]=dfcplot(EEG,EEG.chanlocs,20,10000,500,'Pearson',0.9);
%-----------------------------------------------------------------------------------
%  (CC-BY-4.0) Arthur Valencio
%  Institute of Computing, State University of Campinas
%  Research, Innovation and Dissemination Center for Neuromathematics (RIDC NeuroMat)
%  FAPESP fellowship #2018/0900-8. RIDC NeuroMat also supported by FAPESP #2013/07699-0

%preamble
data=eegdata.data';
timeidx(:,1)=1:step:length(data(:,1));
timeidx(:,2)=timeidx(:,1)+windowlen;
if timeidx(end,2)<=length(eegdata.times)
    times(:,1)=eegdata.times(timeidx(:,1))/1000;
    times(:,2)=eegdata.times(timeidx(:,2))/1000;
else
    idx=find(timeidx(:,2)<=length(eegdata.times),1,'last');
    times(:,1)=eegdata.times(timeidx(1:idx,1))/1000;
    times(:,2)=eegdata.times(timeidx(1:idx,2))/1000;
end
nfigs=ceil(length(times(:,1))/9);
ngrid=rem(length(times(:,1)),9);
k=0;

%calls the DFC fuinctions and calculates the links
[links,lags]=calclinks(data,shift,windowlen,step,dfctype,threshold);

%generate plots
if ~or(strcmp(dfctype,'info'),strcmp(dfctype,'Info'))
    abslinks=abs(links);
    for i=1:nfigs
       figure
       if i==nfigs %last set of (up to) 16 plots
            for j=1:ngrid
                k=k+1;
                try
                    subplot(3,3,j);
                    ds1(k)=topolink(abslinks(k,:,:),chanlocs);
                    title(strcat('t=',num2str(times(k,1),'%.0f'),'-',num2str(times(k,2),'%.0f'),'s'));
                catch
                end
            end
       else %other sets of 16 plots
            for j=1:9
                k=k+1;
                subplot(3,3,j);
                ds1(k)=topolink(abslinks(k,:,:),chanlocs);
                title(strcat('t=',num2str(times(k,1),'%.0f'),'-',num2str(times(k,2),'%.0f'),'s'));
            end
       end
    end
else
    for type=1:6 %different figures for each information measure
        for i=1:nfigs
           figure
           if i==nfigs %last set of (up to) 16 plots
                for j=1:ngrid
                    k=k+1;
                    subplot(3,3,j);
                    ds1(k)=topolink(links(type,k,:,:),chanlocs);
                    title(strcat('t=',num2str(times(k,1),'%.0f'),'-',num2str(times(k,2),'%.0f','s')));
                end
           else
                for j=1:16 %other sets of 16 plots
                    k=k+1;
                    subplot(3,3,j);
                    ds1(k)=topolink(links(type,k,:,:),chanlocs);
                    title(strcat('t=',num2str(times(k,1),'%.0f'),'-',num2str(times(k,2),'%.0f'),'s'));
                end
           end
        end
        %gives a title to figure indicating which information measure it refers to
        if type==1
            sgtitle('CaMI X\rightarrow Y');
        elseif type==2
            sgtitle('CaMI Y\rightarrow X');
        elseif type==3
            sgtitle('Mutual Information');
        elseif type==4
            sgtitle('Transfer Entropy X\rightarrow Y');
        elseif type==5
            sgtitle('Transfer Entropy Y\rightarrow X');
        elseif type==6
            sgtitle('Directionality Index');
        end 
    end
end

end

function [links,lags]=calclinks(data,shift,windowlen,step,dfctype,threshold)
%Calls DFC functions and selects the links above threshold

if or(strcmp(dfctype,'info'),strcmp(dfctype,'Info'))
    %using information-theoretical measures
    %prompt user for required extra variables
    ini_part=input('Initial partitions:');
    L=input('Symbolic length L:');
    units=input('Units:','s');
    %calculates the DFC
    [cami_xy,cami_yx,mutual_info,diridx,te_xy,te_yx] = ...
    dfccami(data,shift,windowlen,step,ini_part,L,units);
    %gets the links
    nintervals=length(cami_xy{1,2}.value);
    nchans=length(data(1,:));
    links(1:6,1:nintervals,1:nchans,1:nchans)=NaN;
    lags(1:6,1:nintervals,1:nchans,1:nchans)=NaN;
    for t=1:nintervals
        for i=1:nchans
            for j=1:nchans
                if i~=j
                  if abs(cami_xy{i,j}.value)>=threshold(1)
                      links(1,t,i,j)=cami_xy{i,j}.value;
                      lags(1,t,i,j)=cami_xy{i,j}.lag;
                  end
                  if abs(cami_yx{i,j}.value)>=threshold(1)
                      links(1,t,i,j)=cami_yx{i,j}.value;
                      lags(1,t,i,j)=cami_yx{i,j}.lag;
                  end
                  if abs(mutual_info{i,j}.value)>=threshold(2)
                      links(1,t,i,j)=mutual_info{i,j}.value;
                      lags(1,t,i,j)=mutual_info{i,j}.lag;
                  end
                  if abs(te_xy{i,j}.value)>=threshold(3)
                      links(1,t,i,j)=te_xy{i,j}.value;
                      lags(1,t,i,j)=mutual_info{i,j}.lag;
                  end
                  if abs(te_yx{i,j}.value)>=threshold(3)
                      links(1,t,i,j)=te_yx{i,j}.value;
                      lags(1,t,i,j)=te_yx{i,j}.lag;
                  end
                  if abs(diridx{i,j}.value)>=threshold(4)
                      links(1,t,i,j)=diridx{i,j}.value;
                      lags(1,t,i,j)=diridx{i,j}.lag;
                  end
                end
            end
        end
    end    
else%using correlation measures
    %calculate DFC
    [corrs,delays] = dfc(data,shift,windowlen,step,dfctype);
    %gets the links
    nintervals=length(corrs(1,1,:));
    nchans=length(data(1,:));
    links(1:nintervals,1:nchans,1:nchans)=NaN;
    lags(1:nintervals,1:nchans,1:nchans)=NaN;
    for t=1:nintervals
       for i=1:nchans
           for j=1:nchans
               if i~=j
                   if abs(corrs(i,j,t))>=threshold
                        links(t,i,j)=corrs(i,j,t);
                        lags(t,i,j)=delays(i,j,t);
                   end
               end
           end
       end
    end
end
end

function [ds]=topolink(links,chanlocs)
%Calls Praneeth Namburi's "Connected Topoplot" script to draw the heads
%Requires EEGLAB

%Creates required structure files
nchans=length(links(1,1,:));
k=0;
for i=1:nchans
    for j=i+1:nchans
        if i~=j
            if ~isnan(links(1,i,j))
                k=k+1;
                ds.chanPairs(k,1)=i;
                ds.chanPairs(k,2)=j;
                ds.connectStrength(k)=links(1,i,j);            
            end
        end
    end
end

%Calls Namburi's draing function
topoplot_connect(ds,chanlocs);

end