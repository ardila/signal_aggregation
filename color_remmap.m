
function remmapped_tensor = color_remmap(M, G)
%Some code to match the output tensor as closely as possible to the
%calculated color dists using a color remapping as described by a matrix,
%ColorM (input_color#, output_color#) which contains weights to linear
%remap each channel

%For now, just calculate according to the first frame
tic
print 'Calculating color distributions'
target_color_dists = color_dists(G.I);

remmapped_tensor = zeros(size(M));
for t = 1:size(M, 4)
    for i = 1:size(M,3)
        remmapped_tensor(:,:,i,t) = ...
            histeq(M(:,:,i,t), target_color_dists(i,:));
    end
end
        