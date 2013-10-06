
%Initialize parameters of generation object, G
G = init(); 


%Create a movie by applying filters to sections, and colorizing according
%to pitch class and spectral roughness.
M = SpectralTensorCrunch(G);

%Remap the color by matching the color distributions in specified by some
%file in G
remmapped_M = color_remmap(M, G);





