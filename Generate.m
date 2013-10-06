%% dash

%Initialize parameters of generation object, G
G = init(); 

%% crunch

%Create a movie by applying filters to sections
M = SpectralTensorCrunch(G);

%Remap the color by matching the color distributions in specified by some
%file in G
remmapped_M = color_remmap(M, G);





