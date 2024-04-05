function score = scoringFunction(truth, reco)


    if any(size(truth) ~= [256 256]) 
        error('The shape of the given ground truth is not 256 x 256!');
    end
    if any(size(reco) ~= [256 256])
        score = 0;
        warning('The shape of the reconstruction is not 256 x 256, so 0 points was given!');
        return;
    end

    truth_c = zeros(size(truth));
    truth_c(truth == 2) = 1;
    reco_c = zeros(size(reco));
    reco_c(reco == 2) = 1;
    
    score_c = KTCssim(reco_c, truth_c);

    truth_d = zeros(size(truth));
    truth_d(truth == 1) = 1;
    reco_d = zeros(size(reco));
    reco_d(reco == 1) = 1;
    
    score_d = KTCssim(reco_d, truth_d);
    

    score = 0.5*(score_c + score_d);

    
   
end





