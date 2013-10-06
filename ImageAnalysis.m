% Get the histogram for each color channel in the image
I = double(imread('Final_Art_2_mb.jpg'));
channels = reshape(I, size(I,3), size(I,1)*size(I,2));
for i = 1:3
    subplot(3,1,i)
    hist(channels(i,:),100)
end