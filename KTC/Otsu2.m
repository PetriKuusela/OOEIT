function [level,x] = Otsu2(image,nvals,~)
% OTSU automatic thresholding for three classes of pixels

[histogramCounts,x] = hist(image(:),nvals);

total = sum(histogramCounts); % total number of pixels in the image

top = 256;

maximum = 0.0;
muT = dot(1:256, histogramCounts(1:end))/sum(histogramCounts);
for ii = 1:top
    for jj = 1:ii
        w1 = sum(histogramCounts(1:jj));
        w2 = sum(histogramCounts(jj+1:ii));
        w3 = sum(histogramCounts(ii+1:end));
        mu1 = (dot(1:jj, histogramCounts(1:jj)))/w1;
        mu2 = (dot(jj+1:ii, histogramCounts(jj+1:ii)))/w2;
        mu3 = (dot(ii+1:256, histogramCounts(ii+1:end)))/w3;
        if w1 > 0 && w2 > 0 && w3 > 0
            val = w1*(mu1-muT)^2 + w2*(mu2-muT)^2 + w3*(mu3-muT)^2;
            if ( val >= maximum )
                level = [jj ii];
                maximum = val;
            end
        end
    end
end
%figure(figno), clf, hist(image(:),256), hold on;
%plot(x(level(1))*ones(2,1),[0,max(histogramCounts)],'LineWidth',2,'Color','r')
%plot(x(level(2))*ones(2,1),[0,max(histogramCounts)],'LineWidth',2,'Color','r')
%title('histogram of image pixels')
%set(gcf,'Units','normalized','OuterPosition',[0.3 0.2 0.3 0.4])
end