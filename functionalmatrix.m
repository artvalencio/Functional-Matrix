function [matrixcorrs,matrixlags] = functionalmatrix(data,shift,windowlen,corrtype)
%FUNCTIONALMATRIX Retrieves the functional matrix from EEG data
%-----------------------------------------------------------------------------------
%Inputs
%- data:        the EEG data, given as a matrix, with rows as time-series and
%               vectors as channels
%- shift:       the maximum shift allowed for shifted correlation analysis, in
%               time-series points units
%- windowlen:   the window length for the correlation analysis (square
%               window)
%- corrtype:    correlation type ('pearson', 'kendall' or 'spearman')
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
    n_frames=length(data(:,1));
    matrixcorrs(1:n_chans,1:n_chans,1:n_frames)=NaN;
    matrixlags(1:n_chans,1:n_chans,1:n_frames)=NaN;
    
    %caluclating correlations and lags
    for i=1:n_chans
        for j=i:n_chans
            [matrixcorrs(i,j,:),matrixlags(i,j,:)]=slidingcorr(data(:,i),data(:,j),shift,windowlen,corrtype);
        end
    end
    
    %complementing triangular matrix without re-doing calcs
    for i=1:n_chans
       for j=1:i
          matrixcorrs(i,j,:)=-matrixcorrs(j,i,:);
          matrixlags(i,j,:)=matrixlags(j,i,:);
       end
    end
end

function [corrval,lag] = slidingcorr(x,y,shift,windowlen,corrtype)
    tslen=length(x);
    %calculates correlations:
    for shift2=-shift:shift
        if shift2>=0
            for i=shift:tslen-windowlen-shift
                r=corr(x(i:i+windowlen),y(i+shift2:i+shift2+windowlen),corrtype);
                correlation(i,shift2+shift+1)=r(1,2);
            end
        else
            for i=shift:tslen-windowlen-shift
                r=corr(x(i+shift2:i+shift2+windowlen),y(i:i+windowlen),corrtype);
                correlation(i,shift2+shift+1)=r(1,2);
            end
        end
    end
    %calculates peak of correlations and respective lags:
    doublecorr=abs(correlation);
    corrlen=length(correlation(:,1));
    corrval(1:corrlen)=NaN;
    lag(a:corrlen)=NaN;
    for i=1:corrlen
        [~,lag(i)]=getpeak(doublecorr(i,:),-shift:shift);
        corrval(i)=correlation(i,lag(i));
    end   
end

function [peak,lag]=getpeak(x,idx)
    peak=-inf;
    for i=2:length(x)
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