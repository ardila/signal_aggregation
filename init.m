function G = init()
    
    %% Image extraction
        

    %% Color reweighting
        G.I = double(imread('Final_Art_2_mb.jpg'));
        G.nbins = 100