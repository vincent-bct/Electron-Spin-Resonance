% datasmooth  Moving average smoothing and differentiation 
%
%   yy = datasmooth(y,m);
%   yy = datasmooth(y,m,'binom')
%   yy = datasmooth(y,m,'flat')
%   yy = datasmooth(y,m,'savgol')
%   yy = datasmooth(y,m,'savgol',p)
%   yy = datasmooth(y,m,'savgol',p,dif)
%
%   Computes the 2*m+1 - point moving average of y.
%   If y is a matrix, datasmooth operates along columns.
%
%   yy(i) is a weighted average of y(i-m:i+m). y is
%   enlarged at start and end by its start and end values,
%   e.g. y(0) = y(-1) = y(1), y(end+2) = y(end).
%
%   If 'flat' is given, the moving average is unweighted.
%
%   If 'binom' is specified, binomial coefficients are
%   used as weighting factors. This is the default method.
%
%   If 'savgol' is specified, a least-squares smoothing using 
%   the Savitzky-Golay polynomial filter of order p is computed.
%   It least-square fits p-order polynomials to 2*m+1 wide frames.
%   If dif>0, a derivative of y is computed at the same time. E.g.
%   if dif=3, y is denoised and its third derivative is returned.
%
%   If the smoothing filter is omitted, 'binom' is the default.
%   If p is omitted, 2 is the default.
%   If dif is omitted, 0 is the default.
