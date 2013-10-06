function M=SpectralTensorCrunch(G) 
%%% Subroutine takes structure G with fields:
% -- fnm:   the name of the audio file to be loaded
% -- Fmn:   the minimum frequency to be computed
% -- Nbnds: the number of frequency bins (must be greater than 2)
% -- Nphs:  the number of phase offsets

%%% A set of basis vectors is calculated. Each basis vector is a sin
%%% function at a different frequency and with a different phase offset
%%% T

%%% Outputs a tensor M of the movie frames as a matrix 
%%% [Frequency x Phase x Time]
fnm=G.fnm;
Fmn=G.Fmn;
Nbnds=G.Nbnds;
Nphs=G.Nphs;
plts=G.plot_spectral_tensor;
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
[fltbnk,ERBff,~]=make_erb_cos_filters(Npts,fs,Nbnds-2,Fmn,fs/2);
% generate the tensor of sinusoid operations
A=zeros(Npts,Nbnds,Nphs);
tt=[0:(Npts-1)]/fs;
for jbnd=1:(Nbnds)
    for jphs=1:Nphs
        A(:,jbnd,jphs)=sin(tt*2*pi*ERBff(jbnd)+(2*pi*(jphs-1)/Nphs));
    end
end

% preallocate a frequency vector for computing harmonicity and roughness
spcff=[1:1:2048]/fs;
N8ve=ceil(log2(fs/2/Fmn));
% scrol through our ERB frequency centers
for jf=1:ERBff
    % for each ERB bin scroll through 8ves and find the values in frequency
    % space of all other bins of the same relative pitch
    ocnt=0;
    for j8ve=-N8ve:N8ve
        [~,ndx]=min(abs(spcff-2^(j8ve)*ERBff(jf)));
        if (ndx~=1 && ndx~=length(spcff))
            ocnt=ocnt+1;
            ERB(jf).octtndx(ocnt)=ndx;
        end
    end
    % and for each ERB bin find the bounds in spcff that correspond to the
    % same ERB pixel
    ERB(jf).bnnddx=find(fltbnk(:,jf)~=0);
end
    
% scroll through sections 
for jsc=1:Nsc
    fprintf('Crunching section %d/%d\n',jsc,Nsc)
    tsc=ts((jsc-1)*Npts/2+[1:Npts],1);
    % generate the spectrum of the time series
    scspc=fft(tsc,4096);
    scspc=scspc(1:2048);
    % scroll through our ER
    % must be done seperately for each channel
    tscL=repmat(tsc,[1 Nbnds Nphs]);
    % perform multiplication and sum over time for power
    M(:,:,jsc,1)=squeeze(sum(tscL.*A,1))/Npts;
    % and for 
    %y(:,:,2,jsc)=squeeze(sum(tscL.*A,1))/Npts;
end
    
% select left channel only for output
M=squeeze(y(:,:,1,:));

rntm=toc;
fprintf('CPUtime=%fs\n',rntm)
%% plot

phhs=[0:1/Nphs:(1-1/Nphs)]*2*pi;

if plts==1
clm=[-100 0];
figure(101)
for jsc=1:Nsc
    subplot(1,3,1)
    imagesc(phhs,ERBff/1e3,20*log10(abs(y(:,:,1,jsc))));
    axis xy
    caxis(clm)
    xlabel('Phase (rad)')
    ylabel('Freq (kHz)')
    title(['Left channel: ' num2str(jsc) '/' num2str(Nsc)])
    subplot(1,3,2)
    imagesc(phhs,ERBff/1e3,20*log10(abs(y(:,:,2,jsc))));
    caxis(clm)
    axis xy
    xlabel('Phase (rad)')
    title(['Right channel: ' num2str(jsc) '/' num2str(Nsc)])
    subplot(1,3,3)
    imagesc(phhs,ERBff/1e3,20*log10(abs(y(:,:,2,jsc)./y(:,:,1,jsc))));
    caxis(clm)
    axis xy
    xlabel('Phase (rad)')
    title(['Right channel: ' num2str(jsc) '/' num2str(Nsc)])
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
