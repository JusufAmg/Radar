clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
R = 110; %initial distance of the target (m)
v = -20; %speed of the target (m/s)

%% FMCW Waveform Generation
Max_Range = 200;
Range_Res = 1;
C = 3e8;
% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
Sweep_Time = 5.5 * 2 * Max_Range / C;
Bandwidth = C / (2*Range_Res);
a = Bandwidth/Sweep_Time;
%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps
%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells
% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Sweep_Time,Nr*Nd); %total time for samples
%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal
%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));
t2=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = R+(v*t(i));
	td(i) = 2*r_t(i)/C;
	t2(i)= t(i)-td(i);
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*((fc*t(i))+(a*t(i)*t(i)/2)));
    Rx (i) = cos(2*pi*((fc*t2(i))+(a*t2(i)*t2(i)/2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix = Rx.*Tx;
    
end

%% RANGE MEASUREMENT

 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix=reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
Mix_fft = fft(Mix,Nr);
 % *%TODO* :
% Take the absolute value of FFT output
Mix_fft=abs(Mix_fft./Nr);
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
L = Sweep_Time * Bandwidth;
Mix_fft=Mix_fft(1:L/2+1);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)
Mix_fft_freq = Bandwidth *(0:(L/2))/L;
plot(Mix_fft_freq,Mix_fft) 
 % *%TODO* :
 % plot FFT output
subplot(2,1,2) 
plot(Mix_fft);
axis ([0 200 0 1]);
%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr=10;
Td=8;
% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr=4;
Gd=4;
% *%TODO* :
% offset the threshold by SNR value in dB
offset=6;
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
Noise_Level = zeros(1,1);
Avg_Noise = zeros(1,1);
Threshold = zeros(1,1);
CFAR = zeros(size(RDM));


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
Grid_Size = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
for ii = Tr+Gr+1:Nr-(Tr+Gr)
    for jj = Td+Gd+1:Nd-(Td+Gd)
        for p = ii-(Tr+Gr) : ii+(Tr+Gr)
            for q = jj-(Td+Gd):jj+(Td+Gd)
                if (abs(ii-p)>Gr || abs(jj-q)>Gd)
                    Noise_Level = Noise_Level+ db2pow(RDM(p,q));
                end
            end
        end
		Avg_Noise = Noise_Level/Grid_Size;
		Threshold = offset +  + pow2db(Avg_Noise);
        if (RDM(ii,jj) >Threshold)
            CFAR(ii,jj) = 1;
        end
    end
end
   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR




% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 








% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,CFAR);
colorbar;


 
 