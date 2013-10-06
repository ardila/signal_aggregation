function color_dists = color_dists(I, n_bins)
% Get the histogram for each color channel in the image I and return it as
% color_dists (channel#, bin#)
channels = reshape(I, size(I,3), size(I,1)*size(I,2));
color_dists = zeros(3, n_bins);
for i = 1:3
    color_dists(i,:) = hist(channels(i,:), n_bins);
end

end


