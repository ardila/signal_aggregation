function M=SpectralTensorCrunch(G) 

%G.fnm='Data/BaktunLoop.wav'; G.Fmn=200; G.Nbnds=256; G.Nphs=256; G.plot_spectral_tensor=1; G.Nft=4096;

%%% Subroutine takes structure G with fields:
% -- fnm:   the name of the audio file to be loaded
% -- Fmn:   the minimum frequency to be computed
% -- Nbnds: the number of frequency bins (must be greater than 2)
% -- Nphs:  the number of phase offsets
% -- Nft:  the number of phase offsets

%%% A set of basis vectors is calculated. Each basis vector is a sin
%%% function at a different frequency and with a different phase offset
%%% T

%%% Outputs a tensor M of the movie frames as a matrix 
%%% [Frequency x Phase x Paramater (intesity, harmonicity, roughness) x Time]
fnm=G.fnm;
Fmn=G.Fmn;
Nbnds=G.Nbnds;
Nphs=G.Nphs;
plts=G.plot_spectral_tensor;
Nft=G.Nft;
%% crunch
tic

% load audio file
[ts,fs]=wavread(fnm);

% divide into 50% overlab sections of length corresponding to one full
% cycle of the lowest frequency.  Ensure this is even.
Npts=ceil(1/Fmn*fs);
Npts=Npts+rem(Npts,2);
% find the number of unique sections that we can fit into the times series
% - zeropad so we can get them all
Nsc=ceil((length(ts)-Npts/2)/Npts);
ts=[ts; zeros(Npts/2*(Nsc+1),size(ts,2))];

% generate the ERB filter centers
[~,ERBff,~]=make_erb_cos_filters(Npts,fs,Nbnds-2,Fmn,fs/2);
% generate the tensor of sinusoid operations FOR THE TIME (call this tensor A)
A=zeros(Npts,Nbnds,Nphs);
tt=[0:(Npts-1)]/fs;
fprintf('Crunching basis vectors, for the time\n')
for jbnd=1:(Nbnds)
    for jphs=1:Nphs
        A(:,jbnd,jphs)=sin(tt*2*pi*ERBff(jbnd)+(2*pi*(jphs-1)/Nphs));
    end
end

% generate the tensor of sinusoid operations, NOT FOR THE TIME but for the
% frequency (call these tensors B and C)
spcff=[1:1:Nft/2]/Nft*fs;
B=zeros(Nft/2,Nbnds,Nphs);
C=zeros(Nft/2,Nbnds,Nphs);
% compute the maximum number of 8ves
N8ve=ceil(log2(fs/2/Fmn));
fprintf('Crunching basis vectors, not for the time, but rather for the frequency\n')
for jbnd=1:(Nbnds)
    % find all values in the frequency spectrum that are 8ves of the thing
    % we care about
    ocnt=0;
    octtndx=[];
    for j8ve=-N8ve:N8ve
        [~,ndx]=min(abs(spcff-2^(j8ve)*ERBff(jbnd)));
        if (ndx~=1 && ndx~=length(spcff))
            ocnt=ocnt+1;
            octtndx(ocnt)=ndx;
        end
    end
    % and now compute all values of frequency within each bin
    bcnt=0;
    bnddx=[];
    for jf=1:Nft/2
        [~,bnndx]=min(abs(spcff(jf)-ERBff));
        if bnndx==jbnd;
            bcnt=bcnt+1;
            bnddx(bcnt)=jf;
        end
    end
    % scroll through phases and compute the necessary tensors
    for jphs=1:Nphs
        B(octtndx,jbnd,jphs)=1*exp(1i*(jphs-1)/Nphs*2*pi);
        C(bnddx,jbnd,jphs)=1*exp(1i*(jphs-1)/Nphs*2*pi);
    end
end
     
% scroll through sections 
M=zeros(Nbnds,Nphs,3,Nsc);
for jsc=1:Nsc
    fprintf('Crunching section %d/%d\n',jsc,Nsc)
    %extract the time series section
    % FOR THE TIME
    tsc=ts((jsc-1)*Npts/2+[1:Npts],1);
    % multiply in time domain and sum over time for power
    tsc=repmat(tsc,[1 Nbnds Nphs]);
    M(:,:,1,jsc)=squeeze(sum(tsc.*A,1))/Npts;
    % FOR THE FREQUENCY generate the spectrum of the time series section
    % and retain only positive freqencies
    scspc=fft(tsc(:,1,1),Nft);
    scspc=scspc(1:Nft/2,:,:);
    scspc=repmat(scspc,[1 Nbnds Nphs]);
    % multiply in freq domain and sum over frequencies for harmonicity 
    M(:,:,2,jsc)=20*log10(abs(squeeze(sum(real(scspc.*B),1)./sum(abs(scspc),1)./sum(abs(B),1))));
    % and for roughness
    M(:,:,3,jsc)=20*log10(abs(squeeze(var(real(scspc.*C),0,1)./sum(abs(C),1))));
    %fprintf('size M: %dx%d%dx%d\n',size(M,1),size(M,2),size(M,3),size(M,4))
    save('gotit.mat','G','M')
end

rntm=toc;
fprintf('CPUtime=%fs\n',rntm)
%% plot

phhs=[0:1/Nphs:(1-1/Nphs)]*2*pi;

if plts==1
clm=[-100 0];
figure(101)
for jsc=1:Nsc
    % intensity
    subplot(1,3,1)
    imagesc(phhs,ERBff/1e3,20*log10(abs(M(:,:,1,jsc))));
    axis xy
    caxis(clm)
    xlabel('Phase (rad)')
    ylabel('Freq (kHz)')
    title(['Intensity: ' num2str(jsc) '/' num2str(Nsc)])
    % harmonicity
    subplot(1,3,2)
    imagesc(phhs,ERBff/1e3,(M(:,:,2,jsc)));
    axis xy
    caxis(clm)
    xlabel('Phase (rad)')
    title(['Harmonicity: ' num2str(jsc) '/' num2str(Nsc)])
    % roughness
    subplot(1,3,3)
    imagesc(phhs,ERBff/1e3,(M(:,:,3,jsc)));
    axis xy
    caxis(clm)
    xlabel('Phase (rad)')
    title(['Roughness: ' num2str(jsc) '/' num2str(Nsc)])
    % colorbar
    hc=colorbar;
    set(hc,'position',[0.92 0.2 0.02 0.6])
    drawnow
    pause(0.01)
end

% figure(102)
% for jbnd=1:(Nbnds)
%     for jphs=1:Nphs
%         subplot('position',[0.05+(jphs-1)*0.9/Nphs 0.05+(jbnd-1)*0.9/(Nbnds+2) 0.9/Nphs 0.9/(Nbnds+2)])
%         plot(tt,A(:,jbnd,jphs));
%         set(gca,'xtick',[],'ytick',[],'ylim',[-1 1],'xlim',[0 max(tt)])
%     end
% end
end
