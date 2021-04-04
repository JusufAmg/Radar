Target Generation and Detection 
In this project, we perform the following tasks:
- Configure the FMCW waveform based on the system requirements.
- Define the range and velocity of target and simulate its displacement.
- For the same simulation loop process the transmit and receive signal to determine the beat signal
- Perform Range FFT on the received signal to determine the Range
- Towards the end, perform the CFAR processing on the output of 2nd FFT to display the target.

![Project_radar_target_generation_and_detection](https://user-images.githubusercontent.com/30601726/113507507-9932ce80-9553-11eb-9557-e76ba560ed02.png)

### Radar System Requirements
|Frequency | 77 Ghz|
:----:|:----:
|Range Resolution | 1 m|
|Max Range | 200 m|
|Max velocity | 70 m/s|
|Velocity resolution| m m/s|

### Target velocity and range
velocity = -20 m/s <br /> 
range = 110 m

### Range FFT (1st FFT)
#### FFT Operation
- Implement the 1D FFT on the Mixed Signal
- Reshape the vector into Nr*Nd array.
- Run the FFT on the beat signal along the range bins dimension (Nr)
- Normalize the FFT output with length, L = Bsweep * Tchirp.
- Take the absolute value of that output.
- Keep one half of the signal
- Plot the output
- There should be a peak at the initial position of the target. In our case at 110 m.

![untitled3](https://user-images.githubusercontent.com/30601726/113507608-3aba2000-9554-11eb-98e4-8b578cd522fd.jpg)



###  2nd FFT for Range Doppler Map
![untitled2](https://user-images.githubusercontent.com/30601726/113507579-1cecbb00-9554-11eb-9634-f8693472a66f.jpg)

### 2D CFAR
```
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
Noise_Level = zeros(size(RDM));
Avg_Noise = zeros(size(RDM));
Threshold = zeros(size(RDM));
CFAR = zeros(size(RDM)); 
```
- Slide the cell under test across the complete matrix. Make sure the CUT has margin for Training and Guard cells from the edges.
- For every iteration sum the signal level within all the training cells. To sum convert the value from logarithmic to linear using db2pow function.
- Average the summed values for all of the training cells used. After averaging convert it back to logarithmic using pow2db.
- Further add the offset to it to determine the threshold.
- Next, compare the signal under CUT against this threshold.
- If the CUT level > threshold assign it a value of 1, else equate it to 0.
- The process above will generate a thresholded block, which is smaller than the Range Doppler Map as the CUTs cannot be located at the edges of the matrix due to the presence of Target and Guard cells. Hence, those cells will not be thresholded. To keep the map size same as it was before CFAR, equate all the non-thresholded cells to 0.
```
Grid_Size = (2*Tr + 2*Gr + 1)*(2*Td + 2*Gd + 1);
for ii = Tr+Gr+1:Nr/2-(Tr+Gr)
    for jj = Td+Gd+1:Nd-(Td+Gd)
        for p = ii-(Tr+Gr) : ii+(Tr+Gr)
            for q = jj-(Td+Gd):jj+(Td+Gd)
                if (abs(ii-p)>Gr || abs(jj-q)>Gd)
                    Noise_Level(ii,jj) = Noise_Level(ii,jj)+ db2pow(RDM(p,q));
                end
            end
        end
		Avg_Noise(ii,jj) = Noise_Level(ii,jj)/Grid_Size;
		Threshold(ii,jj) = offset + pow2db(Avg_Noise(ii,jj));
        if (RDM(ii,jj) >Threshold(ii,jj))
            CFAR(ii,jj) = 1;
        end
    end
end
```
![untitled](https://user-images.githubusercontent.com/30601726/113507542-e7e06880-9553-11eb-8d62-4f19ce0bf407.jpg)

