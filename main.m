function main
%% Read input image
image = double(imread('lena.bmp'));

%% Quantization factor
q_factor = 50;

%% Prepare payload
[col row] = size(image);
payload_maximum=randi([0,1],1,col*row/64);
payload_size = 100;

payload=[payload_maximum(1:payload_size)];

%% Embed
[stream_AC, stream_DC, stream_watermarked_DC, watermarked_qunatized_DCT, original_quantized_DCT,w_reconstructed,reconstructed,failed_flag]=jpegrwdc(image,q_factor,payload);
if failed_flag == 1
    disp('Could not embed the payload')
    pause
end

%% Check recovery
[payload_recovered, recovered_quantized_DCT] = recover_jpegrwdc(watermarked_qunatized_DCT,q_factor);
if not(isequal(recovered_quantized_DCT,original_quantized_DCT))
    disp('Original quantized DCT coefficient could not get recovered')
    if not(isequal(payload_recovered,payload))
        disp('Payload could not get recovered')
    end
    disp('Program will pause')
    pause
end

% file size
disp(['Original JPEG File size is ' num2str(cell_size([stream_AC stream_DC])) ' bits']) 
disp(['Watermarked JPEG File size is ' num2str(cell_size([stream_AC stream_watermarked_DC])) ' bits']) 
disp([num2str(length(payload_recovered)) ' bits are embedded'])
%% Display results
%watermarked image
figure(1)
imshow(uint8(w_reconstructed))

%original JPEG image
figure(2)
imshow(uint8(reconstructed))
pause
end
function [cell_length]=cell_size(A)
cell_length=length(cell2mat(A));
end
