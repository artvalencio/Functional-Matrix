function [matrixcorrs,matrixlags] = functionalmatrix(data,shift,windowlen)
%FUNCTIONALMATRIX Retrieves the functional matrix from EEG data
%-----------------------------------------------------------------------------------
%Inputs
%- data:        the EEG data, given as a matrix, with rows as time-series and
%               vectors as channels
%- shift:       the maximum shift allowed for shifted correlation analysis, in
%               time-series points units
%- windowlen:   the window length for the correlation analysis (square
%               window)
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
    
    n_chans=length(data(1,:));
    matrixcorrs(1:n_chans)=NaN;
    matrixlags(1:n_chans)=NaN;
    for i=1:n_chans
        for j=i:n_chans
            [matrixcorrs(i,j),matrixlags(i,j)]=slidingcorr(data(:,i),data(:,j),shift,windowlen);
        end
    end
    for i=1:n_chans
       for j=1:i
          matrixcorrs(i,j)=-matrixcorrs(j,i);
          matrixlags(i,j)=matrixlags(j,i);
       end
    end
end

function [corrval,lag] = slidingcorr(x,y,shift,windowlen)
    tslen=length(x);
    %calculates correlations:
    for shift2=-shift:shift
        if shift2>=0
            for i=1:tslen-windowlen-shift2
                correlation(i,shift2+shift+1)=corrcoef(x(i:i+windowlen),y(i+shift2:i+shift2+windowlen),'spearman');
            end
        else
            for i=1:tslen-windowlen-shift2
                correlation(i,shift2+shift+1)=corrcoef(x(i+shift2:i+shift2+windowlen),y(i:i+windowlen),'spearman');
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
       if x(i)>peak
            if and(x(i)>x(i-1),x(i)>x(i+1))
                peak=x(i);
                lag=idx(i);
            end
       end
    end
end