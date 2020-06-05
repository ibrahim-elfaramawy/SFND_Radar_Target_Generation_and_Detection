# SFND_Radar_Target_Generation_and_Detection


# TASK 1. FMCW Waveform Design

Based on the radar specification the Bandwidth(B), chirp time (Tchirp) and slope of the chirp (Slope) are calculated as follows:


```
%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

%dres = c / (2*B)
B = c / (2*dres);

%Tchirp = 5.5 * (2*Rmax) /c , this calculation based on the maximum range
%and 5.5 times of the round trip
Tchirp = 5.5 * (2*Rmax) / c;

%Slope
Slope = B / Tchirp;
```


The slope calculated is 2.0455e13


# TASK 2. Simulation Loop

Based on the constant velocity model and given that the initial target position and velocity are Rinit and Vinit respectively, the simulation loop is implemented as follows:


```
%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target r_t(i) for constant velocity.
    r_t(i) = Rinit + Vinit*t(i);
    %Range Estimation equation : Range = (c/2) * taw . where taw is the
    %trip time delay (t_d)
    t_d(i) = (2*r_t(i))/c; 
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i) + Slope * t(i)^2 /2));
    Rx(i) = cos(2*pi*(fc*(t(i) - t_d(i)) + Slope*(t(i)-t_d(i))^2 /2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i) * Rx(i);
    
end
```



# TASK 3. Range FFT (1st FFT)

The range FFT is implemented based on the beat/Mixed frequency (Mix) after being reshaped to a Nr*Nd array as follows:


```
%% RANGE MEASUREMENT

 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
fb = reshape(Mix, [Nr,Nd]); %Beat frequency reshaped
 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
fft_fb = fft(fb,Nr);
fft_fb_normalized = fft_fb/Nr;

 % *%TODO* :
% Take the absolute value of FFT output
fft_output = abs(fft_fb_normalized);
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
fft_oneside = fft_output(1:Nr/2+1);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 
plot(fft_oneside);
axis ([0 200 0 0.5]);
```


The plotted output shows the peak at 111 meters (the initial target position was 110 so the peak matches the initial position).

<img src="images/matlab_1_image.PNG" />


# TASK 4. 2D CFAR

The 2D CFAR is implemented on the output of the 2D FFT operation (Range Doppler Map RDM) as follows:



*   Determine the number of the training and guard cells for each dimension (range and doppler)
*   Determine the value of the offset that will be added later to the threshold
*   Slide the cell under test (CUT) across the complete RDM matrix through a nested loop to make sure the CUT has margin for Training and Guard cells from the edges.
*   Calculate the noise level for every iteration using the sum of the converted signal level value ( from logarithmic to linear using db2pow function) within all the training cells. 
*   Average the summed values (noise level sum / total number of training cells used) and then convert the result back to logarithmic using pow2db.
*   Further add the offset to the average noise level to determine the threshold.
*   Compare the signal under CUT against this threshold.
*   If the CUT level > threshold assign it the max value of RDM , else equate it to 0.
*   Loop through the whole RDM range and doppler dimensions and assign any element that is not included in the previous loop to zero so that the few cells at the edges that are not thresholded are suppressed.

The code implementation for the previous is done as follows:


```
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
offset=10;
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);

% *%TODO* :
% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing CFAR

% Total number of Grid cells
Grid_Size = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);

% Total number of Guard cells
Guard_Size = (2*Gr+1)*(2*Gd+1)-1;

% Total number of training cells 
Training_Size = Grid_Size - Guard_Size - 1;
```


The main loop:


```
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
for i = Tr+Gr+1:(Nr/2)-(Gr+Tr)
    for j = Td+Gd+1:Nd-(Gd+Td)
        noise_level = zeros(1,1);
        for p = i-(Tr+Gr):i+Tr+Gr
            for q = j-(Td+Gd):j+Td+Gd
                if (abs(i-p)>Gr || abs(j-q)>Gd)
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
            end
        end
        threshold = pow2db(noise_level/Training_Size);
        threshold = threshold + offset;
        
        %Measure the signal of the CUT and compare the value against the
        %threshold
        CUT = RDM(i,j);
        if (CUT < threshold)
            RDM(i,j) = 0;
        else
            RDM(i,j) = max(RDM(:));
        end
    end
end
```


The loop to suppress the few cells at the edges that are not thresholded in the main loop:


```
% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
for i =1:Nr/2
    for j=1:Nd
        if (i > (Tr+Gr) && i<((Nr/2)-(Gr+Tr)) && j > (Td+Gd) && j < (Nd-(Gd+Td)))
            continue
        end
        RDM(i,j) = 0;
    end
end
```


Normalize the output and plot it:


```
% Normalize the RDM matrix to hold values between 0 and 1
RDM_normalized = RDM/max(RDM(:));
 
% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RDM_normalized);
colorbar;
```


The output plot:

<img src="images/matlab_2_image.PNG" />

