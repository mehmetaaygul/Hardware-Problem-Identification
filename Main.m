% This code generates a dataset for the paper called "Identification of Distorted RF Components via Deep Multi-Task Learning"

clear all
close all
clc

samples = 1000; %Number of samples

N = 1e3; % Number of symbols
M = 4; % Modulation order
symb = qammod(floor(rand(1,N)*M),M); % Ideal symbols
OSR = 4; % Over sampling rate

output_Filter = zeros(samples,1);
output_Mixer = zeros(samples,1);
output_Phase_shifter = zeros(samples,1);
output_LO = zeros(samples,1);

input_train = zeros(samples,N*OSR);

%%%Probability of distortion
pr_1=0;
pr_2=1;

for ii=1:samples
    ii
    %% IQ Gain
    %Filter
    pr = pr_1 + (pr_2-pr_1) .* rand(1,1);

    if pr>0.5
        a = 0; %Threshold for IQ Gain (compnent is non-distorted)
        b = 0.6; %Threshold for IQ Gain (compnent is non-distorted)
        I_gain = a + (b-a) .* rand(1,1);
        Q_gain = a + (b-a) .* rand(1,1);
        Filter =1;
    else
        a = 0.6; %Threshold for IQ Gain (compnent is distorted)
        b = 1; %Threshold for IQ Gain (compnent is distorted)
        I_gain = a + (b-a) .* rand(1,1);
        Q_gain = a + (b-a) .* rand(1,1);
        Filter =0;
    end

    %% Phase Noise and FO
    %Local oscillator

    pr = pr_1 + (pr_2-pr_1) .* rand(1,1);

    if pr>0.5
        %PN
        a=0; %Threshold for phanse noise (compnent is non-distorted)
        b=60; %Threshold for phanse noise (compnent is distorted)
        p_n = a + (b-a) .* rand(1,1);

        %FO
        a = 0; %Threshold for phanse noise (compnent is non-distorted)
        b = 60; %Threshold for phanse noise (compnent is non-distorted)
        FO_inp = a + (b-a) .* rand(1,1);
        LO = 1;
    else
        a = 60; %Threshold for phanse noise (compnent is distorted)
        b = 90; %Threshold for phanse noise (compnent is distorted)
        p_n = a + (b-a) .* rand(1,1);

        %FO
        a = 60; %Threshold for frequency offset
        b = 100; %Threshold for frequency offset (compnent is distorted)
        FO_inp = a + (b-a) .* rand(1,1);
        LO = 0;
    end

    %% IQ Offset
    % Mixer

    pr = pr_1 + (pr_2-pr_1) .* rand(1,1);

    if pr>0.5
        a = -0.3; %Threshold for IQ offset (compnent is non-distorted)
        b = 0.3;  %Threshold for IQ offset (compnent is non-distorted)
        I_o = a + (b-a) .* rand(1,1);
        Q_o = a + (b-a) .* rand(1,1);
        Mixer = 1;
    else
        a = -0.5; %Threshold for IQ offset (compnent is distorted)
        b = 0.5; %Threshold for IQ offset (compnent is distorted)
        I_o = a + (b-a) .* rand(1,1);
        Q_o = a + (b-a) .* rand(1,1);

        while abs(I_o)<0.3 || abs(Q_o)<0.3
            I_o = a + (b-a) .* rand(1,1);
            Q_o = a + (b-a) .* rand(1,1);
        end
        Mixer =0;
    end
    %% Quadrature Offset
    % Phase shifter
    pr = pr_1 + (pr_2-pr_1) .* rand(1,1);

    if pr>0.5
        a = 0; %Threshold for quadrature offset (compnent is non-distorted)
        b = 0; %Threshold for quadrature offset (compnent is non-distorted)
        qu_o = a + (b-a) .* rand(1,1);
        Phase_shifter = 1;
    else
        a = 60; %Threshold for quadrature offset (compnent is distorted)
        b = 90; %Threshold for quadrature offset (compnent is distorted)
        qu_o = a + (b-a) .* rand(1,1);
        Phase_shifter = 0;
    end

    %% Input Dataset Creation
    N = length(symb);

    noise = randn(1,N*OSR);
    Wn = 1/(length(symb));
    [B,A] = butter(1, Wn);
    dummy = filter(B,A,[1 zeros([1,N])]);
    norm_fac = sqrt(sum(abs(dummy).^2));
    Pn = filter(B,A,noise);
    Pn = Pn/norm_fac;
    Phase_noise = exp(j*Pn*p_n*(pi/180));
    Fs = OSR*N;

    TX_Frame_UpSampled = upsample(symb, OSR);

    n = (0:length(TX_Frame_UpSampled)-1); % Discrete time samples
    symb_fo = TX_Frame_UpSampled./(exp(1i*2*pi*FO_inp*n/Fs)); %Frequency offset
    symb_phase_n = symb_fo.*Phase_noise; %Phase noise
    symb_qu_o = real(symb_phase_n)+(imag(symb_phase_n)*(sind(qu_o)+j*cosd(qu_o))); %Quadrature offset
    symb_i_q = I_gain*real(symb_qu_o) + j*imag(symb_qu_o)*Q_gain ; %I/Q gain
    input(ii,:) = real(symb_i_q) + j*imag(symb_i_q)+I_o+j*Q_o; %I/Q offset

    %% Output Dataset Creation

    output_Filter(ii,:) = Filter;
    output_Mixer(ii,:) = Mixer;
    output_Phase_shifter(ii,:) = Phase_shifter;
    output_LO(ii,:) = LO;
end


