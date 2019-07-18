function [matrixcorrs,matrixlags] = functionalmatrix(data,shift,windowlen,step,corrtype)
%FUNCTIONALMATRIX Retrieves the functional matrix from EEG data
%-----------------------------------------------------------------------------------
%Inputs
%- data:        the EEG data, given as a matrix, with rows as time-series and
%               vectors as channels
%- shift:       the maximum shift allowed for shifted correlation analysis, in
%               time-series points units
%- windowlen:   the window length for the correlation analysis (square
%               window)
%- step:        step given between calculation of correlation in one window
%               and on the other
%- corrtype:    correlation type ('Pearson', 'Kendall' or 'Spearman')
%               (requires Matlab Statistics Toolbox. default: 'Pearson')
%-----------------------------------------------------------------------------------
%Outputs
%- matrixcorrs: Functional matrix of correlations
%- matrixlags:  Matrix of lags to which peak correlations were found
%-----------------------------------------------------------------------------------
%Usage example
%[matrixcorrs,matrixlags] = functionalmatrix(data,20,250)
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
    if and(license('test','statistics_toolbox'),corrtype~=0)
        for i=1:n_chans
            for j=i:n_chans
                disp(strcat(num2str(i),',',num2str(j)))
                [matrixcorrs(i,j,:),matrixlags(i,j,:)]=slidingcorr2(data(:,i),data(:,j),shift,windowlen,step,corrtype);
            end
        end
    else 
        for i=1:n_chans
            for j=i:n_chans
                disp(strcat(num2str(i),',',num2str(j)))
                [matrixcorrs(i,j,:),matrixlags(i,j,:)]=slidingcorr(data(:,i),data(:,j),shift,windowlen,step);
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
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                r=corrcoef(x(i:i+windowlen),y(i+shift2:i+shift2+windowlen));
                correlation(k,shift2+shift+1)=r(1,2);                
            end
        else
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                r=corrcoef(x(i+shift2:i+shift2+windowlen),y(i:i+windowlen));
                correlation(k,shift2+shift+1)=r(1,2);
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

function [peak,lag]=getpeak(x,idx)
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