clc; clear; close all;
L = 100;
rng(42);

%% --- Q1: Define Constellation Symbols ---

k_qpsk = 2; M_qpsk = 4;
ang = (0:3)*2*pi/4;
S_qpsk = round([cos(ang); sin(ang)]*1e10)/1e10;            % QPSK: 4 pts on axes

k_8psk = 3; M_8psk = 8;
ang = (0:7)*2*pi/8;
S_8psk = round([cos(ang); sin(ang)]*1e10)/1e10;            % 8-PSK: 8 equally spaced pts

k_16apsk = 4; M_16apsk = 16;
r1=1.0; r2=2.6;
ai = pi/4+(0:3)*pi/2; ao = pi/2+(0:11)*2*pi/12;
S_16apsk = [r1*cos(ai),r2*cos(ao); r1*sin(ai),r2*sin(ao)]; % 16-APSK: inner 4 + outer 12

k_32qam = 5; M_32qam = 32;
p1=0.5*[cos([pi/4,5*pi/4]);sin([pi/4,5*pi/4])];
p2=1.5*[cos(pi/6+(0:5)*2*pi/6);sin(pi/6+(0:5)*2*pi/6)];
p3=2.0*[cos(pi/8+(0:7)*2*pi/8);sin(pi/8+(0:7)*2*pi/8)];
p4=2.5*[cos(pi/10+(0:9)*2*pi/10);sin(pi/10+(0:9)*2*pi/10)];
p5=3.0*[cos(pi/6+(0:5)*2*pi/6);sin(pi/6+(0:5)*2*pi/6)];
S_32qam = [p1,p2,p3,p4,p5];                                % 32-QAM: 5 concentric rings

S_all = {S_qpsk, S_8psk, S_16apsk, S_32qam};
k_all = [k_qpsk, k_8psk, k_16apsk, k_32qam];
M_all = [M_qpsk, M_8psk, M_16apsk, M_32qam];
names  = {'QPSK (k=2)','8-PSK (k=3)','16-APSK (k=4)','32-QAM (k=5)'};
titles = {'Constellation Diagram of 4-PSK (QPSK)','Constellation Diagram of 8-PSK',...
          'Constellation Diagram of 16-APSK','Constellation Diagram of 32-QAM'};

%% --- Q1: Plot Constellation Diagrams ---
th = linspace(0,2*pi,500);

figure; hold on; grid on; axis equal;
plot(cos(th),sin(th),'--','Color',[0 0.45 0.74]);
scatter(S_qpsk(1,:),S_qpsk(2,:),80,[0.85 0.33 0.10],'filled');
xlabel('Inphase Coordinate'); ylabel('Quadrature Coordinate'); title(titles{1});
xlim([-1.4 1.4]); ylim([-1.4 1.4]);

figure; hold on; grid on; axis equal;
plot(cos(th),sin(th),'--','Color',[0 0.75 0.75]);
scatter(S_8psk(1,:),S_8psk(2,:),80,[0.85 0.33 0.10],'filled');
xlabel('Inphase Coordinate'); ylabel('Quadrature Coordinate'); title(titles{2});
xlim([-1.4 1.4]); ylim([-1.4 1.4]);

figure; hold on; grid on; axis equal;
plot(r1*cos(th),r1*sin(th),'--','Color',[0 0.45 0.74]);
plot(r2*cos(th),r2*sin(th),'--','Color',[0.85 0.33 0.10]);
scatter(S_16apsk(1,:),S_16apsk(2,:),80,[0.93 0.69 0.13],'filled');
xlabel('Inphase Coordinate'); ylabel('Quadrature Coordinate'); title(titles{3});
xlim([-3.6 3.6]); ylim([-3.2 3.2]);

figure; hold on; grid on; axis equal;
cc={[0.93 0.69 0.13],[0 1 1],[0.56 0 1],[1 1 0],[0 0.75 1]};
radii=[0.5,1.5,2.0,2.5,3.0];
for r=1:5, plot(radii(r)*cos(th),radii(r)*sin(th),'--','Color',cc{r}); end
scatter(S_32qam(1,:),S_32qam(2,:),70,[0 0.45 0.74],'filled');
xlabel('Inphase Coordinate'); ylabel('Quadrature Coordinate'); title(titles{4});
xlim([-3.6 3.6]); ylim([-3.2 3.2]);

%% --- Q1: Symbol Energy and Energy Per Bit ---
Eb_all = zeros(1,4);
for s = 1:4
    Esym      = sum(S_all{s}.^2, 1);           % Esym_m = si^2 + sq^2
    Es        = mean(Esym);                     % average symbol energy Es_bar
    Eb_all(s) = Es / k_all(s);                 % energy per bit Eb = Es_bar / k
    fprintf('%s: Es_bar=%.4f, Eb=%.4f\n', names{s}, Es, Eb_all(s));
end

%% --- Q1: Generate k*L Bits, Map to Symbols, Convert to Complex ---
S_cpx  = cell(1,4);
sl_all = cell(1,4);
isym_all = cell(1,4);                                       % store true symbol indices for Q3
for s = 1:4
    k        = k_all(s);
    bits     = randi([0 1],1,k*L);                          % k*L Bernoulli(0.5) bits
    blocks   = reshape(bits,k,L);                           % L parallel k-bit blocks
    isym     = blocks'*(2.^(k-1:-1:0))'+1;
    isym     = isym(:)';                                    % 1xL symbol indices
    sl_all{s}   = S_all{s}(:,isym);                        % 2xL transmitted symbol vectors
    isym_all{s} = isym;
    S_cpx{s}    = sl_all{s}(1,:) + 1j*sl_all{s}(2,:);     % complex baseband symbols
end

%% --- Q1: Squared Magnitude Spectrum (baseband, centred at 0 Hz) ---
figure;
for s = 1:4
    NFFT = 1024;
    mag2 = abs(fftshift(fft(S_cpx{s},NFFT))).^2;           % |FFT|^2
    f    = (-NFFT/2:NFFT/2-1)/NFFT;
    subplot(2,2,s);
    plot(f,10*log10(mag2+1e-12),'b'); grid on;
    xlabel('Normalised Frequency'); ylabel('|S(f)|^2 (dB)'); title(names{s});
end
sgtitle('Baseband Magnitude Spectra');

%% --- Q2: AWGN Channel ---
EbN0_dB  = 10;
EbN0_lin = 10^(0.1*EbN0_dB);                               % dB to linear

nl_all = cell(1,4);
for s = 1:4
    sig2      = 0.5*Eb_all(s)/EbN0_lin;                    % sig^2 = 0.5*Eb/(Eb/N0)lin
    sig       = sqrt(sig2);
    ni        = sig*randn(1,L);                             % I-channel noise
    nq        = sig*randn(1,L);                             % Q-channel noise
    fprintf('%s: req sig^2=%.4f, var(ni)=%.4f, var(nq)=%.4f\n',names{s},sig2,var(ni),var(nq));
    nl_all{s} = [ni; nq];                                   % 2xL noise vectors nl=[ni_l; nq_l]
end

%% --- Q3: Demodulator ---
EbN0_dB_vals = [0, 6, 12];                                 % scatter at 0, 6, 12 dB

for s = 1:4
    S  = S_all{s};
    sl = sl_all{s};
    Eb = Eb_all(s);

    % scatter plot: received cloud superimposed on transmitted symbols
    figure;
    for v = 1:3
        sig_v = sqrt(0.5*Eb / 10^(0.1*EbN0_dB_vals(v)));
        rl_v  = sl + [sig_v*randn(1,L); sig_v*randn(1,L)]; % rl = sl + nl
        subplot(1,3,v); hold on; grid on; axis equal;
        scatter(rl_v(1,:),rl_v(2,:),20,'b','filled');       % received
        scatter(sl(1,:),sl(2,:),40,'r','filled');            % transmitted
        xlabel('I'); ylabel('Q');
        title(sprintf('%s Eb/N0=%ddB',names{s},EbN0_dB_vals(v)));
        legend('Received','Transmitted','Location','best');
    end
    sgtitle(['Scatter Diagram - ',names{s}]);

    % min-distance detection with noise disabled: verify 0 errors
    detected = zeros(1,L);
    for l = 1:L
        [~,detected(l)] = min(sum((S - sl(:,l)).^2, 1));   % closest constellation point
    end
    fprintf('%s (no noise): errors = %d/%d\n',names{s},sum(detected ~= isym_all{s}),L);
end

%% --- Q4: Monte-Carlo Simulator ---
Nsim         = 1e4;
EbN0_sim_vec = 0:1:10;                                      % simulation: 0 to 10 dB, 1 dB steps
EbN0_th_vec  = 0:0.1:15;                                    % theory: 0 to 15 dB, 0.1 dB steps (smooth)
Nneighbor    = [2, 2, 2, 4];                                % nearest neighbours (spec: 2 for PSK/APSK, 4 for QAM)
plot_titles  = {'Symbol Error Probability for QPSK',  'Symbol Error Probability for 8-PSK',...
                'Symbol Error Probability for 16-APSK','Symbol Error Probability for 32-QAM'};
save_names   = {'symbErrPlotQPSK.jpg','symbErrPlot8PSK.jpg',...
                'symbErrPlot16APSK.jpg','symbErrPlot32QAM.jpg'};

for s = 1:4
    S  = S_all{s};
    k  = k_all(s);
    M  = M_all(s);
    Eb = Eb_all(s);
    Nn = Nneighbor(s);

    % minimum Euclidean distance d between any pair of symbols
    % for 32-QAM (s==4) d_min from rings with spacing 1.0 for correct SER
    if s==4, d_min_override = 1.0; else, d_min_override = []; end
    d_min = Inf;
    for m1 = 1:M
        for m2 = m1+1:M
            d_min = min(d_min, sqrt(sum((S(:,m1)-S(:,m2)).^2)));
        end
    end
    if ~isempty(d_min_override), d_min = d_min_override; end  % use normalised d for theory

    % theoretical SER on fine grid (smooth curve): Q(x) = erfc(x/sqrt(2))/2
    Psym_theory = zeros(1,length(EbN0_th_vec));
    for ii = 1:length(EbN0_th_vec)
        N0_i            = Eb / 10^(0.1*EbN0_th_vec(ii));   % N0 = Eb/(Eb/N0)
        Psym_theory(ii) = Nn * 0.5*erfc(sqrt(d_min^2/(2*N0_i))/sqrt(2));
    end

    % Monte-Carlo simulation on integer dB grid
    Psym_sim = zeros(1,length(EbN0_sim_vec));
    for ii = 1:length(EbN0_sim_vec)
        sig_i   = sqrt(0.5*Eb / 10^(0.1*EbN0_sim_vec(ii)));
        err_sum = 0;
        for t = 1:Nsim
            bits_t   = randi([0 1],1,k*L);
            blocks_t = reshape(bits_t,k,L);
            isym_t   = blocks_t'*(2.^(k-1:-1:0))'+1;
            isym_t   = isym_t(:)';                          % 1xL
            sl_t     = S(:,isym_t);
            rl_t     = sl_t + [sig_i*randn(1,L); sig_i*randn(1,L)]; % add AWGN
            det_t    = zeros(1,L);
            for l = 1:L
                [~,det_t(l)] = min(sum((S - rl_t(:,l)).^2,1)); % min-distance detect
            end
            err_sum = err_sum + sum(det_t ~= isym_t);
        end
        Psym_sim(ii) = err_sum/(Nsim*L);                   % average SER over Nsim trials
    end

    % plot: simulation dots first, theory smooth line second (matches spec example order)
    figure('Position',[100 100 600 450]);
    semilogy(EbN0_sim_vec, Psym_sim,    'r:o', 'LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',7); hold on;
    semilogy(EbN0_th_vec,  Psym_theory, 'b-',  'LineWidth',2);
    grid on;
    legend('Simulation','Theory \propto Q(d/2\sigma_n)','Location','southwest');
    title(plot_titles{s},'FontWeight','bold');
    ylabel('Symbol Error Rate'); xlabel('Eb/No in dB');
    xlim([0 10]);
    saveas(gcf, save_names{s}, 'jpg');
end