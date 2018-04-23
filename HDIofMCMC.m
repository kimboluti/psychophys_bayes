function HDIlim = HDIofMCMC(sampleVec , credMass )
%
% Computes highest density interval from a sample of representative values,
% estimated as shortest credible interval.
% Arguments:
% sampleVec
% is a vector of representative values from a probability distribution.
% credMass
% is a scalar between 0 and 1, indicating the mass within the credible
% interval that is to be estimated.
% Value:
%   HDIlim is a vector containing the limits of the HDI

% sortedPts = sort( sampleVec );
% ciIdxInc = floor( credMass .* length( sortedPts ) );
% nCIs = length( sortedPts ) - ciIdxInc; %identifies the number of possible HDI there can be
% 
% for i = 1:nCIs 
%    ciWidth(i) = sortedPts(i + ciIdxInc ) - sortedPts( i );
% end;
% HDImin = sortedPts( min( ciWidth ) );
% HDImax = sortedPts( HDImin + ciIdxInc );
% HDIlim = [HDImin HDImax];


%computes highest density interval from a sample of representative values, 
%adapted from Kruschke (2011)
%estimated as shortest credible interval
%arguments
% sampleVec
%   is a vector or matrix of representative values from a probability
%   distribution
%credMass
%   is a scalar between 0 and 1 indicateing the mass within the credible
%   interval that is to be estimated
% HDI is a vector with containing the limits of the HDI

sortedPts = sort(sampleVec(:)); % sort the the sampled points, allowing for a matrix or vector with :
ciIdxInc = floor(credMass.*length(sortedPts)); % find the index increment for the HDI
nCIs = length(sortedPts) - ciIdxInc; %how many intervals are possible
i = 1:nCIs; % the possible lower interval estimate index
j = i+ciIdxInc; % the possible higher interval estimates;
ciWidth(i) = sortedPts(j)-sortedPts(i); %calculate the candidate widths;
[theMin, theMinI] = min(ciWidth); %find the index of the minimum

HDIlim(1) = sortedPts(i(theMinI));%the lower interval estimate
HDIlim(2) = sortedPts(j(theMinI));%the upper interval estimate