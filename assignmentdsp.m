training_files_zero_male = dir('C:\Users\Abd\Desktop\DSPassingment\Train\44100M\*.wav');
testing_files_zero_male = dir('C:\Users\Abd\Desktop\DSPassingment\Test\44100M\*.wav');
training_files_zero_female = dir('C:\Users\Abd\Desktop\DSPassingment\Train\44100F\*.wav');
testing_files_zero_female = dir('C:\Users\Abd\Desktop\DSPassingment\Test\44100F\*.wav');

% Total number of male and female test samples
total_male = length(testing_files_zero_male);
sum_signal_male = 0;
total_female = length(testing_files_zero_female);
sum_signal_female = 0;

% Determine the maximum length of the audio files
max_length = 0;
for i = 1:length(training_files_zero_male)
    file_path = strcat(training_files_zero_male(i).folder, '\', training_files_zero_male(i).name);
    y = audioread(file_path);
    max_length = max(max_length, length(y));
end

avg_psd_male = zeros(floor(max_length/2)+1, 1);
% Initialize sum_signal_male and sum_signal_female with zeros
sum_signal_male = zeros(max_length, 1);

% ------------ Training -----------------------------

% read the 'zero male' training files and calculate the energy of them.
data_zero_male_energy = [];
data_zero_male_ZCR=[];
for i = 1:length(training_files_zero_male)
file_path = strcat(training_files_zero_male(i).folder,'\',training_files_zero_male(i).name);
[y,fs] = audioread(file_path);

% Zero-pad the signal if necessary
    if length(y) < max_length
        y = [y; zeros(max_length - length(y), 1)];
    end

ZCR_zero_male_train1 = mean(abs(diff(sign(y(1:floor(end/3))))))./2;
ZCR_zero_male_train2 = mean(abs(diff(sign(y(floor(end/3):floor (end*2/3))))))./2;
ZCR_zero_male_train3 = mean(abs(diff(sign(y(floor(end*2/3):end)))))./2;

energy_zero_male=sum(y .^2);
ZCR_zero_male = [ZCR_zero_male_train1 ZCR_zero_male_train2 ZCR_zero_male_train3];
sum_signal_male = sum_signal_male + y;

data_zero_male_energy = [data_zero_male_energy energy_zero_male];
data_zero_male_ZCR=[data_zero_male_ZCR ; ZCR_zero_male];


    y = [y; zeros(max_length - length(y), 1)]; % Zero-padding to max_length
    Y = fft(y);
    P2 = abs(Y/length(y)).^2;
    P1 = P2(1:length(y)/2+1);
    avg_psd_male = avg_psd_male + P1;
end

energy_zero_male=mean(data_zero_male_energy);
fprintf('The avarege  energy of male saying zero is \n');
disp(energy_zero_male);

ZCR_zero_male= mean( data_zero_male_ZCR);
fprintf('ZCR of Male Saying Zero is \n');
disp(ZCR_zero_male);
avg_signal_male = sum_signal_male / length(training_files_zero_male);


avg_psd_male = avg_psd_male / length(training_files_zero_male);


%-----------------------------------------------------------------------

for i = 1:length(training_files_zero_female)
    file_path = strcat(training_files_zero_female(i).folder, '\', training_files_zero_female(i).name);
    y = audioread(file_path);
    max_length = max(max_length, length(y));
end
sum_signal_female = zeros(max_length, 1);
avg_psd_female = zeros(floor(max_length/2)+1, 1);
% read the 'zero female' training files and calculate the energy of them.
data_zero_female_energy = [];
data_zero_female_ZCR=[];
for i = 1:length(training_files_zero_female)
file_path = strcat(training_files_zero_female(i).folder,'\',training_files_zero_female(i).name);
[y,fs] = audioread(file_path);

% Zero-pad the signal if necessary
    if length(y) < max_length
        y = [y; zeros(max_length - length(y), 1)];
    end
    
ZCR_zero_female_train1 = mean(abs(diff(sign(y(1:floor(end/3))))))./2;
ZCR_zero_female_train2 = mean(abs(diff(sign(y(floor(end/3):floor (end*2/3))))))./2;
ZCR_zero_female_train3 = mean(abs(diff(sign(y(floor(end*2/3):end)))))./2;

energy_zero_female=sum(y .^2);
ZCR_zero_female = [ZCR_zero_female_train1 ZCR_zero_female_train2 ZCR_zero_female_train3];

data_zero_female_energy = [data_zero_female_energy energy_zero_female];
data_zero_female_ZCR=[data_zero_female_ZCR ; ZCR_zero_female];
sum_signal_female = sum_signal_female + y;

    y = [y; zeros(max_length - length(y), 1)]; % Zero-padding to max_length
    Y = fft(y);
    P2 = abs(Y/length(y)).^2;
    P1 = P2(1:length(y)/2+1);
    avg_psd_female = avg_psd_female + P1;


end

energy_zero_female=mean(data_zero_female_energy);
fprintf('The energy of Female saying zero is \n');
disp(energy_zero_female);

ZCR_zero_female= mean(data_zero_female_ZCR);
fprintf('ZCR of female Saying Zero is \n');
disp(ZCR_zero_female);

avg_signal_female = sum_signal_female / length(training_files_zero_female);


avg_psd_female = avg_psd_female / length(training_files_zero_female);



% ------------ Evaluation -----------------------------
fprintf('----------Evaluating According to Energy----------\n');

% Initialize counters for Energy
tp_energy_male = 0;
fp_energy_male = 0;
tp_energy_female = 0;
fp_energy_female = 0;

% read the 'male' tesing files and calculate the energy of them.
for i = 1:length(testing_files_zero_male)
file_path = strcat(testing_files_zero_male(i).folder,'\',testing_files_zero_male(i).name);
[y,fs] = audioread(file_path);

y_energy  = sum(y.^2);
    %>|
    % Plot in Time Domain
    figure;
    subplot(2,1,1); % Subplot 1 for time domain
    t = linspace(0, length(y)/fs, length(y)); % Time vector
    plot(t, y);
    title(['Test file [Zero Male] #' num2str(i) ' in Time Domain']);
    xlabel('Time (seconds)');
    ylabel('Amplitude');

    % Plot in Frequency Domain
    Y = fft(y); % Fourier Transform
    P2 = abs(Y/length(y));
    P1 = P2(1:floor(length(y)/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(length(y)/2))/length(y); % Frequency vector
    subplot(2,1,2); % Subplot 2 for frequency domain
    plot(f, P1);
    title(['Test file [Zero Male] #' num2str(i) ' in Frequency Domain']);
    xlabel('Frequency (Hz)');
    ylabel('|P1(f)|');
       %^
 % test if the energy of this file is closer to YES or NO average energies
    if(abs(y_energy-energy_zero_male) < abs(y_energy-energy_zero_female)) 
        fprintf('Test file [Zero Male] #%d classified as Male saying Zero ,E=%d\n',i,y_energy);
        tp_energy_male = tp_energy_male + 1;
    else
        fprintf('Test file [Zero Male] #%d classified as Female Saying Zero E=%d\n',i,y_energy);
        
    end
end
%---------------------------------------------------
% read the 'female' tesing files and calculate the energy of them.
for i = 1:length(testing_files_zero_female)
file_path = strcat(testing_files_zero_female(i).folder,'\',testing_files_zero_female(i).name);
[y,fs] = audioread(file_path);

y_energy  = sum(y.^2);

%>|
    % Plot in Time Domain
    figure;
    subplot(2,1,1); % Subplot 1 for time domain
    t = linspace(0, length(y)/fs, length(y)); % Time vector
    plot(t, y);
    title(['Test file [Zero Female] #' num2str(i) ' in Time Domain']);
    xlabel('Time (seconds)');
    ylabel('Amplitude');

    % Plot in Frequency Domain
    Y = fft(y); % Fourier Transform
    P2 = abs(Y/length(y));
    P1 = P2(1:floor(length(y)/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(length(y)/2))/length(y); % Frequency vector
    subplot(2,1,2); % Subplot 2 for frequency domain
    plot(f, P1);
    title(['Test file [Zero Female] #' num2str(i) ' in Frequency Domain']);
    xlabel('Frequency (Hz)');
    ylabel('|P1(f)|');
       %^

    if(abs(y_energy-energy_zero_male) < abs(y_energy-energy_zero_female))
        fprintf('Test file [Zero Female] #%d classified as Female Saying Zero ,E=%d\n',i,y_energy);
        tp_energy_female = tp_energy_female + 1;

        
    else
        fprintf('Test file [Zero Female] #%d classified as Male Saying Zero ,E=%d\n',i,y_energy);
        

    end
end
% Calculate accuracy
accuracy_energy_male = tp_energy_male / total_male;
accuracy_energy_female = tp_energy_female / total_female;
% Display the accuracies
fprintf('Accuracy for male classification using Energy: %f\n', accuracy_energy_male);
fprintf('Accuracy for female classification using Energy: %f\n', accuracy_energy_female);

fprintf('\n----------Evaluating According to ZCR----------\n');
% Initialize counters for ZCR
tp_ZCR_male = 0;
fp_ZCR_male = 0;
tp_ZCR_female = 0;
fp_ZCR_female = 0;

% read the 'male' tesing files and calculate the ZCR of them.
for i = 1:length(testing_files_zero_male)
file_path = strcat(testing_files_zero_male(i).folder,'\',testing_files_zero_male(i).name);
[y,fs] = audioread(file_path);

%divide the signal into 3 parts and calculate the ZCR for each part
ZCR_male1 = mean(abs(diff(sign(y(1:floor(end/3))))))./2;
ZCR_male2 = mean(abs(diff(sign(y(floor(end/3):floor (end*2/3))))))./2;
ZCR_male3 = mean(abs(diff(sign(y(floor(end*2/3):end)))))./2;

y_ZCR = [ZCR_male1 ZCR_male2 ZCR_male3];
 % test if the energy of this file is closer to YES or NO average energies
    if(pdist([y_ZCR;ZCR_zero_male],'cosine') < pdist([y_ZCR;ZCR_zero_female],'cosine'))
        fprintf('Test file [Zero Male] #%d classified as Male saying Zero \n',i);
        tp_ZCR_male=tp_ZCR_male+1;
    else
        fprintf('Test file [Zero Male] #%d classified as Female Saying Zero\n',i);
        
        
    end
end

% read the 'female' tesing files and calculate the ZCR of them.

for i = 1:length(testing_files_zero_female)
file_path = strcat(testing_files_zero_female(i).folder,'\',testing_files_zero_female(i).name);
[y,fs] = audioread(file_path);

%divide the signal into 3 parts and calculate the ZCR for each part
ZCR_female1 = mean(abs(diff(sign(y(1:floor(end/3))))))./2;
ZCR_female2 = mean(abs(diff(sign(y(floor(end/3):floor (end*2/3))))))./2;
ZCR_female3 = mean(abs(diff(sign(y(floor(end*2/3):end)))))./2;

y_ZCR = [ZCR_female1 ZCR_female2 ZCR_female3];
 % test if the energy of this file is closer to YES or NO average energies
    if(pdist([y_ZCR;ZCR_zero_male],'cosine') < pdist([y_ZCR;ZCR_zero_female],'cosine'))
        fprintf('Test file [Zero female] #%d classified as female saying Zero \n',i);
        tp_ZCR_female=tp_ZCR_female+1;
    else
        fprintf('Test file [Zero female] #%d classified as male Saying Zero\n',i);
        
    end
end
accuracy_ZCR_male = tp_ZCR_male / total_male;
accuracy_ZCR_female = tp_ZCR_female / total_female;
fprintf('Accuracy for male classification using ZCR: %f\n', accuracy_ZCR_male);
fprintf('Accuracy for female classification using ZCR: %f\n', accuracy_ZCR_female);


%--------------------------------------------------------------------------------
fprintf('\n----------Evaluating According to Correlation----------\n');
% Initialize counters for ZCR
tp_corr_male = 0;
fp_corr_male = 0;
tp_corr_female = 0;
fp_corr_female = 0;

% read the 'male' tesing files and calculate the correlation of them.
for i = 1:length(testing_files_zero_male)
    file_path = strcat(testing_files_zero_male(i).folder, '\', testing_files_zero_male(i).name);
    [y, fs] = audioread(file_path);

    % Make y the same length as avg_signal_male
    if length(y) < length(avg_signal_male)
        y_padded_male = [y; zeros(length(avg_signal_male) - length(y), 1)];
    else
        y_padded_male = y(1:length(avg_signal_male));
    end

    % Compute correlation with male average
    corr_with_male = corrcoef(y_padded_male, avg_signal_male);
    corr_value_male = corr_with_male(1,2);

    % Make y the same length as avg_signal_female
    if length(y) < length(avg_signal_female)
        y_padded_female = [y; zeros(length(avg_signal_female) - length(y), 1)];
    else
        y_padded_female = y(1:length(avg_signal_female));
    end

    % Compute correlation with female average
    corr_with_female = corrcoef(y_padded_female, avg_signal_female);
    corr_value_female = corr_with_female(1,2);

    % Classify based on which correlation is higher
    if corr_value_male > corr_value_female
        fprintf('Test file [Zero Male] #%d classified as Male saying Zero, Correlation=%f\n', i, corr_value_male);
        tp_corr_male = tp_corr_male + 1;
    else
        fprintf('Test file [Zero Male] #%d classified as Female Saying Zero, Correlation=%f\n', i, corr_value_female);
        
    end
end


% read the 'female' tesing files and calculate the correlation of them.
for i = 1:length(testing_files_zero_female)
    file_path = strcat(testing_files_zero_female(i).folder, '\', testing_files_zero_female(i).name);
    [y, fs] = audioread(file_path);

    % Make y the same length as avg_signal_male
    if length(y) < length(avg_signal_male)
        y_padded_male = [y; zeros(length(avg_signal_male) - length(y), 1)];
    else
        y_padded_male = y(1:length(avg_signal_male));
    end

    % Compute correlation with male average
    corr_with_male = corrcoef(y_padded_male, avg_signal_male);
    corr_value_male = corr_with_male(1,2);

    % Make y the same length as avg_signal_female
    if length(y) < length(avg_signal_female)
        y_padded_female = [y; zeros(length(avg_signal_female) - length(y), 1)];
    else
        y_padded_female = y(1:length(avg_signal_female));
    end

    % Compute correlation with female average
    corr_with_female = corrcoef(y_padded_female, avg_signal_female);
    corr_value_female = corr_with_female(1,2);

    % Classify based on which correlation is higher
    if corr_value_male < corr_value_female
        fprintf('Test file [Zero Male] #%d classified as Male saying Zero, Correlation=%f\n', i, corr_value_male);
        tp_corr_female = tp_corr_female + 1;
    else
        fprintf('Test file [Zero Male] #%d classified as Female Saying Zero, Correlation=%f\n', i, corr_value_female);
        
    end
end

accuracy_crr_male = tp_corr_male / total_male;
accuracy_crr_female = tp_corr_female / total_female;
fprintf('Accuracy for male classification using Correlation: %f\n', accuracy_crr_male);
fprintf('Accuracy for female classification using Correlation: %f\n', accuracy_crr_female);

%----------------------------------------------------------------------------------------------------------------

fprintf('\n----------PSD Evaluating----------\n');
% Initialize counters for ZCR
tp_psd_male = 0;
tp_psd_female = 0;


% Analyze PSD for male testing files
for i = 1:length(testing_files_zero_male)
    file_path = strcat(testing_files_zero_male(i).folder, '\', testing_files_zero_male(i).name);
    [y, fs] = audioread(file_path);

    % Pad or truncate y to match the longer of avg_psd_male or avg_psd_female
    max_psd_length = max(length(avg_psd_male), length(avg_psd_female));
    if length(y) < max_psd_length
        y_padded = [y; zeros(max_psd_length - length(y), 1)]; % Zero padding
    else
        y_padded = y(1:max_psd_length); % Truncating
    end

    % Compute the PSD of the padded signal
    Y = fft(y_padded);
    P2 = abs(Y/length(y_padded)).^2;
    P1 = P2(1:floor(length(y_padded)/2)+1); % Ensure P1 is half the FFT result

    % Compare with average PSDs
    % Ensure P1 is the same length as avg_psd_male for comparison
    if length(P1) < length(avg_psd_male)
        P1_padded = [P1; zeros(length(avg_psd_male) - length(P1), 1)]; % Zero padding
    else
        P1_padded = P1(1:length(avg_psd_male)); % Truncating
    end
    distance_to_male = sum((P1_padded - avg_psd_male).^2);

    % Ensure P1 is the same length as avg_psd_female for comparison
    if length(P1) < length(avg_psd_female)
        P1_padded = [P1; zeros(length(avg_psd_female) - length(P1), 1)]; % Zero padding
    else
        P1_padded = P1(1:length(avg_psd_female)); % Truncating
    end
    distance_to_female = sum((P1_padded - avg_psd_female).^2);

    % Classification logic
    if distance_to_male < distance_to_female
        fprintf('Test file [Zero Male] #%d classified as Male\n', i);
        tp_psd_male=tp_psd_male+1;
    else
        fprintf('Test file [Zero Male] #%d classified as Female\n', i);
    end
end

%-------------------------------------------------------------------------------------

for i = 1:length(testing_files_zero_female)
    file_path = strcat(testing_files_zero_female(i).folder, '\', testing_files_zero_female(i).name);
    [y, fs] = audioread(file_path);

    % Pad or truncate y to match the longer of avg_psd_male or avg_psd_female
    max_psd_length = max(length(avg_psd_male), length(avg_psd_female));
    if length(y) < max_psd_length
        y_padded = [y; zeros(max_psd_length - length(y), 1)]; % Zero padding
    else
        y_padded = y(1:max_psd_length); % Truncating
    end

    % Compute the PSD of the padded signal
    Y = fft(y_padded);
    P2 = abs(Y/length(y_padded)).^2;
    P1 = P2(1:floor(length(y_padded)/2)+1); % Ensure P1 is half the FFT result

    % Compare with average PSDs
    % Ensure P1 is the same length as avg_psd_male for comparison
    if length(P1) < length(avg_psd_male)
        P1_padded = [P1; zeros(length(avg_psd_male) - length(P1), 1)]; % Zero padding
    else
        P1_padded = P1(1:length(avg_psd_male)); % Truncating
    end
    distance_to_male = sum((P1_padded - avg_psd_male).^2);

    % Ensure P1 is the same length as avg_psd_female for comparison
    if length(P1) < length(avg_psd_female)
        P1_padded = [P1; zeros(length(avg_psd_female) - length(P1), 1)]; % Zero padding
    else
        P1_padded = P1(1:length(avg_psd_female)); % Truncating
    end
    distance_to_female = sum((P1_padded - avg_psd_female).^2);

    % Classification logic
    if distance_to_male < distance_to_female
        fprintf('Test file [Zero Female] #%d classified as female\n', i);
        tp_psd_female=tp_psd_female+1;
    else
        fprintf('Test file [Zero female] #%d classified as male\n', i);
    end
end
accuracy_psd_male = tp_psd_male / total_male;
accuracy_psd_female = tp_psd_female / total_female;
fprintf('Accuracy for male classification using PSD: %f\n', accuracy_psd_male);
fprintf('Accuracy for female classification using PSD: %f\n', accuracy_psd_female);

%----------------------Merge classifications: Energy ZCR PSD Correlation

fprintf('-----------------Merge Classifications---------------\n');
tp_merge_male = 0;
fp_merge_male = 0;
tp_merge_female = 0;
fp_merge_female = 0;
for i = 1:length(testing_files_zero_male)
    file_path = strcat(testing_files_zero_male(i).folder, '\', testing_files_zero_male(i).name);
    [y, fs] = audioread(file_path);
    y_energy  = sum(y.^2);
    
    ZCR_male1 = mean(abs(diff(sign(y(1:floor(end/3))))))./2;
    ZCR_male2 = mean(abs(diff(sign(y(floor(end/3):floor (end*2/3))))))./2;
    ZCR_male3 = mean(abs(diff(sign(y(floor(end*2/3):end)))))./2;
    y_ZCR = [ZCR_male1 ZCR_male2 ZCR_male3];
    
    
     % Make y the same length as avg_signal_male
    if length(y) < length(avg_signal_male)
        y_padded_male = [y; zeros(length(avg_signal_male) - length(y), 1)];
    else
        y_padded_male = y(1:length(avg_signal_male));
    end

    % Compute correlation with male average
    corr_with_male = corrcoef(y_padded_male, avg_signal_male);
    corr_value_male = corr_with_male(1,2);

    % Make y the same length as avg_signal_female
    if length(y) < length(avg_signal_female)
        y_padded_female = [y; zeros(length(avg_signal_female) - length(y), 1)];
    else
        y_padded_female = y(1:length(avg_signal_female));
    end

    % Compute correlation with female average
    corr_with_female = corrcoef(y_padded_female, avg_signal_female);
    corr_value_female = corr_with_female(1,2);
    
     % Pad or truncate y to match the longer of avg_psd_male or avg_psd_female
    max_psd_length = max(length(avg_psd_male), length(avg_psd_female));
    if length(y) < max_psd_length
        y_padded = [y; zeros(max_psd_length - length(y), 1)]; % Zero padding
    else
        y_padded = y(1:max_psd_length); % Truncating
    end

    % Compute the PSD of the padded signal
    Y = fft(y_padded);
    P2 = abs(Y/length(y_padded)).^2;
    P1 = P2(1:floor(length(y_padded)/2)+1); % Ensure P1 is half the FFT result

    % Compare with average PSDs
    % Ensure P1 is the same length as avg_psd_male for comparison
    if length(P1) < length(avg_psd_male)
        P1_padded = [P1; zeros(length(avg_psd_male) - length(P1), 1)]; % Zero padding
    else
        P1_padded = P1(1:length(avg_psd_male)); % Truncating
    end
    distance_to_male = sum((P1_padded - avg_psd_male).^2);

    % Ensure P1 is the same length as avg_psd_female for comparison
    if length(P1) < length(avg_psd_female)
        P1_padded = [P1; zeros(length(avg_psd_female) - length(P1), 1)]; % Zero padding
    else
        P1_padded = P1(1:length(avg_psd_female)); % Truncating
    end
    distance_to_female = sum((P1_padded - avg_psd_female).^2);
    
    if((abs(y_energy-energy_zero_male) < abs(y_energy-energy_zero_female)) || (pdist([y_ZCR;ZCR_zero_male],'cosine') < pdist([y_ZCR;ZCR_zero_female],'cosine'))||(corr_value_male > corr_value_female) || (distance_to_male < distance_to_female)) 
        fprintf('Test file [Zero Male] #%d classified as Male saying Zero \n',i);
        tp_merge_male=tp_merge_male+1;
    else
        fprintf('Test file [Zero Male] #%d classified as Female Saying Zero\n',i);
        
    end
end


%female testing 
for i = 1:length(testing_files_zero_female)
    file_path = strcat(testing_files_zero_female(i).folder, '\', testing_files_zero_female(i).name);
    [y, fs] = audioread(file_path);
    y_energy  = sum(y.^2);
    
    ZCR_female1 = mean(abs(diff(sign(y(1:floor(end/3))))))./2;
    ZCR_female2 = mean(abs(diff(sign(y(floor(end/3):floor (end*2/3))))))./2;
    ZCR_female3 = mean(abs(diff(sign(y(floor(end*2/3):end)))))./2;
    y_ZCR = [ZCR_female1 ZCR_female2 ZCR_female3];
    
    
     % Make y the same length as avg_signal_male
    if length(y) < length(avg_signal_male)
        y_padded_male = [y; zeros(length(avg_signal_male) - length(y), 1)];
    else
        y_padded_male = y(1:length(avg_signal_male));
    end

    % Compute correlation with male average
    corr_with_male = corrcoef(y_padded_male, avg_signal_male);
    corr_value_male = corr_with_male(1,2);

    % Make y the same length as avg_signal_female
    if length(y) < length(avg_signal_female)
        y_padded_female = [y; zeros(length(avg_signal_female) - length(y), 1)];
    else
        y_padded_female = y(1:length(avg_signal_female));
    end

    % Compute correlation with female average
    corr_with_female = corrcoef(y_padded_female, avg_signal_female);
    corr_value_female = corr_with_female(1,2);
    
     % Pad or truncate y to match the longer of avg_psd_male or avg_psd_female
    max_psd_length = max(length(avg_psd_male), length(avg_psd_female));
    if length(y) < max_psd_length
        y_padded = [y; zeros(max_psd_length - length(y), 1)]; % Zero padding
    else
        y_padded = y(1:max_psd_length); % Truncating
    end

    % Compute the PSD of the padded signal
    Y = fft(y_padded);
    P2 = abs(Y/length(y_padded)).^2;
    P1 = P2(1:floor(length(y_padded)/2)+1); % Ensure P1 is half the FFT result

    % Compare with average PSDs
    % Ensure P1 is the same length as avg_psd_male for comparison
    if length(P1) < length(avg_psd_male)
        P1_padded = [P1; zeros(length(avg_psd_male) - length(P1), 1)]; % Zero padding
    else
        P1_padded = P1(1:length(avg_psd_male)); % Truncating
    end
    distance_to_male = sum((P1_padded - avg_psd_male).^2);

    % Ensure P1 is the same length as avg_psd_female for comparison
    if length(P1) < length(avg_psd_female)
        P1_padded = [P1; zeros(length(avg_psd_female) - length(P1), 1)]; % Zero padding
    else
        P1_padded = P1(1:length(avg_psd_female)); % Truncating
    end
    distance_to_female = sum((P1_padded - avg_psd_female).^2);
    
    if((abs(y_energy-energy_zero_male) < abs(y_energy-energy_zero_female)) || (pdist([y_ZCR;ZCR_zero_male],'cosine') < pdist([y_ZCR;ZCR_zero_female],'cosine'))||(corr_value_male > corr_value_female) || (distance_to_male < distance_to_female)) 
        fprintf('Test file [Zero female] #%d classified as female saying Zero \n',i);
        tp_merge_female=tp_merge_female+1;
    else
        fprintf('Test file [Zero female ] #%d classified as male Saying Zero\n',i);
        
    end
end


accuracy_merge_male = tp_merge_male / total_male;
accuracy_merge_female = tp_merge_female / total_female;
fprintf('Accuracy for male classification using all classifications: %f\n', accuracy_merge_male);
fprintf('Accuracy for female classification using all classifications: %f\n', accuracy_merge_female);

