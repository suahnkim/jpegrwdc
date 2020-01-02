function main
%% Read input image
image_name='lena512';
image = double(imread([image_name '.bmp']));

%% Quantization factor
q_factor = 70;

%% Prepare payload
[col row] = size(image);
payload_maximum=randi([0,1],1,col*row/64);
payload_size = 500;

payload=[payload_maximum(1:payload_size)];

%% Embed
[original_JPEG, watermarked_JPEG,failed_flag]=jpegrwdc_pixel(image,q_factor,payload);
if failed_flag == 1
    disp('Could not embed the payload')
    pause
end
jpeg_write(original_JPEG.JPEG_struct,[image_name '.jpg']);
jpeg_write(watermarked_JPEG.JPEG_struct,[image_name '_watermarked.jpg']);

%% Check recovery
watermarked_image = jpeg_read([image_name '_watermarked.jpg']);
[payload_recovered, recovered_JPEG] = recover_jpegrwdc(watermarked_image);
jpeg_write(recovered_JPEG,[image_name '_recovered.jpg']);

original=double(imread([image_name '.jpg']));
recovered=double(imread([image_name '_recovered.jpg']));
watermarked=double(imread([image_name '_watermarked.jpg']));

if not(isequal(original,recovered))
    disp('Original quantized DCT coefficient could not get recovered')
    disp('Program will pause')
    pause
else
    disp('Original quantized DCT coefficient are recovered')
end
if not(isequal(payload_recovered,payload))
    disp('Payload could not get recovered')
    disp('Program will pause')
    pause
else
    disp('Payload has been recovered')
end

original_file=dir([image_name '.jpg']);
watermarked_file=dir([image_name '_watermarked.jpg']);
% file size
disp(['Original JPEG File size is ' num2str(original_file.bytes) ' bits']) 
disp(['Watermarked JPEG File size is ' num2str(watermarked_file.bytes) ' bits']) 
disp([num2str(length(payload_recovered)) ' bits are embedded'])
%% Display results
%watermarked image
figure(1)
imshow(uint8(watermarked))

%original JPEG image
figure(2)
imshow(uint8(original))
end
% function [cell_length]=cell_size(A)
% cell_length=length(cell2mat(A));
% end
