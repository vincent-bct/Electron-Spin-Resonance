% basecorr  Baseline correction 
%
%   CorrSpectrum = basecorr(Spectrum,Dimension,Order)
%   [CorrSpectrum,BaseLine] = basecorr(...)
%
%   Makes a polynomial baseline correction for the array Spectrum.
%   Dimension gives a list of dimensions to which successive 1D
%   corrections should be applied. If Dimensions=[], then a multi-
%   dimensional (hyper)surface fitting is performed along all
%   dimensions.
%   Order gives a list of the polynomial orders for all dimensions
%   specified. All orders must be between 0 and 4.
%
%   Examples:
%      z = basecorr(z,1,2); % 2nd order along dimension 1
%      z = basecorr(z,[],3); % twodimensional 3rd order surface
%      z = basecorr(z,[1 2],[3 1]); % 3rd order along dim 1 and 1st
%                                   % order along dim 2.
