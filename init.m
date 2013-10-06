function G = init()
    addpath('Utils')
    addpath('Data')
    %% Image extraction
        %%% Subroutine takes structure G with fields:
        % -- fnm:   the name of the audio file to be loaded
        % -- Fmn:   the minimum frequency to be computed
        % -- Nbnds: the number of frequency bins (must be greater than 2)
        % -- Nphs:  the number of phase offsets

        %%% Outputs a tensor M of the movie frames as a matrix 
        %%% [Frequency x Phase x Time]

        G.fnm='BaktunLoop.wav';
        G.Fmn=200;
        G.Nbnds=256;
        G.Nphs=256;
        G.plot_spectral_tensor = 0;
        
    %% Colorization
        %%% After "image extraction" there is one array, I. We want color in
        %%% our final result. That means that for each pixel, we must
        %%% gather two additional pieces of information to supplement the
        %%% value already calculated and allow us to embed each pixel in a 
        %%% 3 dimensional color space.
        %%% In this implementation, we use the total power in a pitch class
        %%% (C1, C2, C3, C4) summed across all octaves as one piece of
        %%% information.
        %%% The second piece of information is spectral roughness. Since
        %%% each pixel corresponds to a frequency bin which contains
        %%% several different frequencies, it is possible to measure the
        %%% power at each of the frequencies within the bin and then
        %%% calculate the variance of these values: a measure of how much
        %%% power varies within that frequency bin.
        %%% Now, each pixel has 3 pieces of information which can be used
        %%% to embed each pixel in a 3 dimensional color space. RGB is
        %%% often used, but we are using LAB because it more closely
        %%% approximated the representation at lower levels of the human
        %%% visual system.
 
    

    %% Color reweighting
        %%%I: image to use as target for color distributions
        %%%nbins: number of bins to use in calculating color distributions
        G.I = double(imread('Final_Art_5_meg.jpg'));
        G.nbins = 100;