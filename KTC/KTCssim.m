function score = KTCssim(truth, reco, r)
%This function calculates a structural similarity index measure based score for the
%reconstruction "reco", comparing it to ground truth "truth". "nc" is the
%number of different classes of pixels, i.e. the number of different values
%the pixels have.
%
%The scoring is based on the structural similarity index measure:
%https://en.wikipedia.org/wiki/Structural_similarity
%The notations are almost similar to the ones used in the Wikipedia article

    if any(size(truth) ~= size(reco))
        error('The ground truth and reconstruction have different sizes!');
    end
    if nargin < 3 || isempty(r)
        r = 80;
    end

    c1 = 1e-4;
    c2 = 9e-4;
    
    %create a gaussian filter kernel
    ws = ceil(2*r);
    [X, Y] = meshgrid(-ws:ws);
    ker = 1./(r*sqrt(2*pi))*exp(-0.5*(X.^2 + Y.^2)./r^2);
    correction = conv2(ones(size(truth)), ker, 'same');

    gt = conv2(truth, ker, 'same')./correction;
    gr = conv2(reco, ker, 'same')./correction;
    
    mu_t2 = gt.^2;
    mu_r2 = gr.^2;
    mu_t_mu_r = gt.*gr;

    sigma_t2 = conv2(truth.^2, ker, 'same')./correction - mu_t2;
    sigma_r2 = conv2(reco.^2, ker, 'same')./correction - mu_r2;
    sigma_tr = conv2(truth.*reco, ker, 'same')./correction - mu_t_mu_r;

    num = (2*mu_t_mu_r + c1).*(2*sigma_tr + c2);
    den = (mu_t2 + mu_r2 + c1).*(sigma_t2 + sigma_r2 + c2);
    ssimimage = num./den;

    score = mean(ssimimage(:));

end



