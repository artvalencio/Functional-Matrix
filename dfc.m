function [matrixcorrs,matrixlags] = dfc(data,shift,windowlen,step,corrtype)
%DFC Retrieves the Dynamic Functional Connectivity matrix from EEG data
%using correlation measures
%-----------------------------------------------------------------------------------
%Inputs
%- data:        the EEG data, given as a matrix, with each column being a
%               EEG channel
%- shift:       the maximum shift allowed for correlation analysis considering 
%               time-delay between locations. Given in time-series points
%               units. Select shift=0 for disregard time-delay effects
%- windowlen:   the interval (window length) on which the correlation is 
%               calculated
%- step:        step (jump) given between calculation of correlation in one window
%               and on the other. Select step=1 for better precision,
%               step>1 for faster computation
%- corrtype:    correlation type ('Pearson', 'Kendall' or 'Spearman')
%               (for other options that not 'Pearson', it requires 
%               Matlab Statistics Toolbox or Octave Statistics Package. 
%               default: 'Pearson')
%
%NOTE: TO USE INFORMATION MEASURES INSTEAD, USE DFCCAMI
%-----------------------------------------------------------------------------------
%Outputs
%- matrixcorrs: Functional matrix of correlations
%- matrixlags:  Matrix of lags to which peak correlations were found
%-----------------------------------------------------------------------------------
%Usage example
%[matrixcorrs,matrixlags] = dfc(data,20,250,1,'Pearson')
%-----------------------------------------------------------------------------------
%  (CC-BY-4.0) Arthur Valencio
%  Institute of Computing, State University of Campinas
%  Research, Innovation and Dissemination Center for Neuromathematics (RIDC NeuroMat)
%  FAPESP fellowship #2018/0900-8. RIDC NeuroMat also supported by FAPESP #2013/07699-0
    
    %initializing
    n_chans=length(data(1,:));
    
    %check if corrtype
    if nargin<5
        corrtype=0;
    end
    
    %calculating correlations and lags
    if ~isOctave
        if and(license('test','statistics_toolbox'),corrtype~=0)
            disp(strcat('Calculating ',corrtype,' correlations'));
            for i=1:n_chans
                for j=i:n_chans
                    disp(strcat(num2str(i),',',num2str(j)))
                    [matrixcorrs(i,j,:),matrixlags(i,j,:)]=slidingcorr2(data(:,i),data(:,j),shift,windowlen,step,corrtype);
                end
            end
        else 
            disp(strcat('Calculating Pearson correlations'));
            for i=1:n_chans
                for j=i:n_chans
                    disp(strcat(num2str(i),',',num2str(j)))
                    [matrixcorrs(i,j,:),matrixlags(i,j,:)]=slidingcorr(data(:,i),data(:,j),shift,windowlen,step);
                end
            end
        end
    else
        if or(corrtype==0,or(corrtype=='Pearson',corrtype=='pearson'))
            disp(strcat('Calculating Pearson correlations'));
            for i=1:n_chans
                for j=i:n_chans
                    disp(strcat(num2str(i),',',num2str(j)))
                    [matrixcorrs(i,j,:),matrixlags(i,j,:)]=slidingcorr(data(:,i),data(:,j),shift,windowlen,step);
                end
            end
        else
            try
                disp(strcat('Calculating ',corrtype,' correlations'));
                for i=1:n_chans
                    for j=i:n_chans
                        disp(strcat(num2str(i),',',num2str(j)))
                        [matrixcorrs(i,j,:),matrixlags(i,j,:)]=slidingcorr3(data(:,i),data(:,j),shift,windowlen,step,corrtype);
                    end
                end
            catch
                try
                    pkg load statistics
                    disp(strcat('Calculating ',corrtype,' correlations'));
                    for i=1:n_chans
                        for j=i:n_chans
                            disp(strcat(num2str(i),',',num2str(j)))
                            [matrixcorrs(i,j,:),matrixlags(i,j,:)]=slidingcorr3(data(:,i),data(:,j),shift,windowlen,step,corrtype);
                        end
                    end
                catch
                    disp(strcat('Calculating Pearson correlations'));
                    for i=1:n_chans
                        for j=i:n_chans
                            disp(strcat(num2str(i),',',num2str(j)))
                            [matrixcorrs(i,j,:),matrixlags(i,j,:)]=slidingcorr(data(:,i),data(:,j),shift,windowlen,step);
                        end
                    end
                end
            end
        end
    end
    %complementing triangular matrix without re-doing calcs
    for i=1:n_chans
       for j=1:i-1
          matrixcorrs(i,j,:)=-matrixcorrs(j,i,:);
          matrixlags(i,j,:)=matrixlags(j,i,:);
       end
    end
end

function [corrval,lag] = slidingcorr(x,y,shift,windowlen,step)
    tslen=length(x);
    %calculates correlations:
    for shift2=-shift:shift
        if shift2>=0
            k=0;
            %disp(tslen-windowlen-shift)
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                r=corrcoef(x(i:i+windowlen),y(i+shift2:i+shift2+windowlen));
                correl(k,shift2+shift+1)=r(1,2);                
            end
        else
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                r=corrcoef(x(i+shift2:i+shift2+windowlen),y(i:i+windowlen));
                correl(k,shift2+shift+1)=r(1,2);
            end
        end
    end
    %calculates peak of correlations and respective lags:
    abscorr=abs(correl);
    corrlen=length(correl(:,1));
    corrval(1:corrlen)=NaN;
    lag(1:corrlen)=NaN;
    for i=1:corrlen
        if ~isempty(find(abscorr(i,:),1)) %if correlation is not totally zero
            [~,lag(i)]=getpeak(abscorr(i,:),-shift:shift);
            corrval(i)=correl(i,lag(i)+shift+1);
        else %else correltion is zero and lag undefined
            lag(i)=NaN;
            corrval(i)=0;
        end
    end   
end

function [corrval,lag] = slidingcorr2(x,y,shift,windowlen,step,corrtype)
    tslen=length(x);
    %calculates correlations:
    for shift2=-shift:shift
        if shift2>=0
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                correlation(k,shift2+shift+1)=corr(x(i:i+windowlen),...
                    y(i+shift2:i+shift2+windowlen),'Type',corrtype,'Rows','Pairwise');
            end
        else
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                correlation(k,shift2+shift+1)=corr(x(i+shift2:i+shift2+windowlen),...
                    y(i:i+windowlen),'Type',corrtype,'Rows','Pairwise');
            end
        end
    end
    %calculates peak of correlations and respective lags:
    abscorr=abs(correlation);
    corrlen=length(correlation(:,1));
    corrval(1:corrlen)=NaN;
    lag(1:corrlen)=NaN;
    for i=1:corrlen
        if ~isempty(find(abscorr(i,:),1)) %if correlation is not totally zero
            [~,lag(i)]=getpeak(abscorr(i,:),-shift:shift);
            corrval(i)=correlation(i,lag(i)+shift+1);
        else %else correltion is zero and lag undefined
            lag(i)=NaN;
            corrval(i)=0;
        end
    end   
end

function [corrval,lag] = slidingcorr3(x,y,shift,windowlen,step,corrtype)
    tslen=length(x);
    %calculates correlations:
    %-case Spearman
    if or(corrtype=='Spearman',corrtype=='spearman')
    for shift2=-shift:shift
        if shift2>=0
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                correlation(k,shift2+shift+1)=spearman(x(i:i+windowlen),...
                    y(i+shift2:i+shift2+windowlen));
            end
        else
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                correlation(k,shift2+shift+1)=spearman(x(i+shift2:i+shift2+windowlen),...
                    y(i:i+windowlen));
            end
        end
    end
    %-case Kendall
    elseif or(corrtype=='kendall',or(corrtype=='Kendall',or(corrtype=='kendal',corrtype=='Kendal')))
        if shift2>=0
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                correlation(k,shift2+shift+1)=kendall(x(i:i+windowlen),...
                    y(i+shift2:i+shift2+windowlen));
            end
        else
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                correlation(k,shift2+shift+1)=kendall(x(i+shift2:i+shift2+windowlen),...
                    y(i:i+windowlen));
            end
        end
    end
    %calculates peak of correlations and respective lags:
    abscorr=abs(correlation);
    corrlen=length(correlation(:,1));
    corrval(1:corrlen)=NaN;
    lag(1:corrlen)=NaN;
    for i=1:corrlen
        if ~isempty(find(abscorr(i,:),1)) %if correlation is not totally zero
            [~,lag(i)]=getpeak(abscorr(i,:),-shift:shift);
            corrval(i)=correlation(i,lag(i)+shift+1);
        else %else correltion is zero and lag undefined
            lag(i)=NaN;
            corrval(i)=0;
        end
    end   
end

function [peak,lag]=getpeak(x,idx)
%find peak of absolute value of correlation within -shift:shift interval
    peak=-inf;
    for i=2:length(x)-1
       if ~isnan(x(i))
            if x(i)>peak
                if and(x(i)>x(i-1),x(i)>x(i+1))
                    peak=x(i);
                    lag=idx(i);
                end
            end
       end
    end
end