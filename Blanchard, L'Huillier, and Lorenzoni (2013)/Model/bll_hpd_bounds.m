function [lower, upper] = bll_hpd_bounds(draws, coverage)
%BLL_HPD_BOUNDS Shortest empirical posterior intervals along last dimension.

arguments
    draws double {mustBeFinite}
    coverage (1,1) double {mustBeGreaterThan(coverage,0),mustBeLessThan(coverage,1)}
end

originalSize = size(draws);
drawCount = originalSize(end);
if drawCount < 2
    error('BLL:InsufficientDraws', 'At least two posterior draws are required.');
end
rowCount = numel(draws)/drawCount;
sortedDraws = sort(reshape(draws, rowCount, drawCount), 2);
intervalSteps = max(1, min(drawCount-1, ceil(coverage*drawCount)-1));
widths = sortedDraws(:,intervalSteps+1:end) ...
    - sortedDraws(:,1:end-intervalSteps);
[~, startIndex] = min(widths, [], 2);
rows = (1:rowCount)';
lowerVector = sortedDraws(sub2ind(size(sortedDraws), rows, startIndex));
upperVector = sortedDraws(sub2ind(size(sortedDraws), rows, startIndex+intervalSteps));
summarySize = originalSize(1:end-1);
lower = reshape(lowerVector, summarySize);
upper = reshape(upperVector, summarySize);
end
