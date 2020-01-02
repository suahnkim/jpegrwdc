% JPEG Reversible Watermarking in quantized DC coefficient
% Suah Kim
% Created: November 2012
% Last Modified: December 2019

function [original_JPEG, watermarked_JPEG, failed_flag]=jpegrwdc_pixel(image,q_factor,payload)


failed_flag=0;
DCTQ_50=[...
    16 11 10 16 24 40 51 61;...
    12 12 14 19 26 58 60 55;...
    14 13 16 24 40 57 69 56;...
    14 17 22 29 51 87 80 62;...
    18 22 37 56 68 109 103 77;...
    24 36 55 64 81 194 113 92;...
    49 64 78 87 103 121 120 101;...
    72 92 95 98 112 100 103 99];

order_zig=[1,2,9,17,10,3,4,11,18,25,33,26,19,12,5,6,13,20,27,34,41,49,42,35,28,21,14,7,8,15,22,29,36,43,50,57,58,51,44,37,30,23,16,24,31,38,45,52,59,60,53,46,39,32,40,47,54,61,62,55,48,56,63,64];
%% Read Image

[height, width]=size(image);
if q_factor <50
    s_factor = 5000/q_factor;
else
    s_factor = 200 - q_factor*2;
end


DCTQ_50 = DCTQ_50*s_factor + 50;
DCTQ_50 = (DCTQ_50/100);


for i = 1:1:8
    for j = 1:1:8
        if DCTQ_50(i,j) <1
            DCTQ_50(i,j) = 1;
        end
    end
end
DCTQ=double(uint8(round(DCTQ_50)));
%% Check Size for JPEG compliance
if mod(width,2)~=0
    width = width-1;
end

if mod(height,2)~=0
    height = height-1;
end

image = imresize(image, [height width]);
%% Initialize Variables
%Load JPEG tables for initialization (JPEG_struct)
original_JPEG=load('JPEG_grey_init_table.mat');
watermarked_JPEG=load('JPEG_grey_init_table.mat');

original_JPEG.JPEG_struct.width= width;
original_JPEG.JPEG_struct.height= height;
original_JPEG.JPEG_struct.quant_tables=mat2cell(DCTQ,8,8);

watermarked_JPEG.JPEG_struct.width= width;
watermarked_JPEG.JPEG_struct.height= height;
watermarked_JPEG.JPEG_struct.quant_tables=mat2cell(DCTQ,8,8);

%% Convert Image To YCBCR
d= size(image);
if d == 3
    img = rgb2ycbcr(image);
else
    img=image(:,:,1);
end

%% Segmenting Image Blocks Of 8x8
k=0;
Quantized_DCT=zeros(size(img,1),size(img,2));
dctmtx_8=dctmtx(8);
for index_i=1:size(img,1)/8
    for index_j=1:size(img,2)/8
        i=1+8*(index_i-1);
        j=1+8*(index_j-1);
        Quantized_DCT(i:i+7,j:j+7)=round(dctmtx_8*(img(i:i+7,j:j+7)-128)*dctmtx_8'./DCTQ);
        reconstructed(i:i+7,j:j+7)=round(dctmtx_8'*(Quantized_DCT(i:i+7,j:j+7).*DCTQ)*dctmtx_8)+128;
    end
end

[sorted_list,zig_zag_coeff]=prediction_list(Quantized_DCT,DCTQ,order_zig,0);

%% Embedding
%Initialize
embedded_payload_size=0;
payload_length_size_binary=ceil(log2(height*width/64));
flag=0;
watermarked_qc=Quantized_DCT;
payload_size=length(payload);
side_information_size=de2bi(payload_size,payload_length_size_binary);
embedded_side_information_size=0;

%Embed even set
for i=1:length(sorted_list)
    if flag==0
        dc_prediction=sorted_list(i,3);
        i1=sorted_list(i,1);
        j1=sorted_list(i,2);
        if dc_prediction == 0
            embedded_payload_size=embedded_payload_size+1;
            watermarked_qc(i1,j1)=watermarked_qc(i1,j1)+payload(embedded_payload_size);
            if embedded_payload_size == payload_size/2
                flag=1;
            end
        elseif dc_prediction == -1
            embedded_payload_size=embedded_payload_size+1;
            watermarked_qc(i1,j1)=watermarked_qc(i1,j1)-payload(embedded_payload_size);
            if embedded_payload_size == payload_size/2
                flag=1;
            end
        elseif dc_prediction > 0
            watermarked_qc(i1,j1)=watermarked_qc(i1,j1)+1;
        elseif dc_prediction < -1
            watermarked_qc(i1,j1)=watermarked_qc(i1,j1)-1;
        end
    end
end
flag=0;

%Embed odd set
[sorted_list_odd,~]=prediction_list(watermarked_qc,DCTQ,order_zig,1);

for i=1:length(sorted_list_odd)
    if flag==0
        dc_prediction=sorted_list_odd(i,3);
        i1=sorted_list_odd(i,1);
        j1=sorted_list_odd(i,2);
        if dc_prediction == 0
            if embedded_side_information_size < payload_length_size_binary
                embedded_side_information_size=embedded_side_information_size+1;
                watermarked_qc(i1,j1)=watermarked_qc(i1,j1)+side_information_size(embedded_side_information_size);
            else
                embedded_payload_size=embedded_payload_size+1;
                watermarked_qc(i1,j1)=watermarked_qc(i1,j1)+payload(embedded_payload_size);
                if embedded_payload_size == payload_size
                    flag=1;
                end
            end
        elseif dc_prediction == -1
            if embedded_side_information_size < payload_length_size_binary
                embedded_side_information_size=embedded_side_information_size+1;
                watermarked_qc(i1,j1)=watermarked_qc(i1,j1)-side_information_size(embedded_side_information_size);
            else
                embedded_payload_size=embedded_payload_size+1;
                watermarked_qc(i1,j1)=watermarked_qc(i1,j1)-payload(embedded_payload_size);
                if embedded_payload_size == payload_size
                    flag=1;
                end
            end
        elseif dc_prediction > 0
            watermarked_qc(i1,j1)=watermarked_qc(i1,j1)+1;
        elseif dc_prediction < -1
            watermarked_qc(i1,j1)=watermarked_qc(i1,j1)-1;
        end
    end
    if  i == length(sorted_list_odd) && flag == 0
        flag = 2;
    end
end

if flag ==2
    disp('Payload could not be embeddd')
    failed_flag=1;
end

counter=0;


watermarked_DC=zeros(size(Quantized_DCT,1)/64*size(Quantized_DCT,2),1);
for index_i=1:size(img,1)/8
    for index_j=1:size(img,2)/8
        i2=1+8*(index_i-1);
        j2=1+8*(index_j-1);
        counter=counter+1;
        watermarked_DC(counter)=watermarked_qc(i2,j2);
        w_reconstructed(i2:i2+7,j2:j2+7)=round(dctmtx_8'*(watermarked_qc(i2:i2+7,j2:j2+7).*DCTQ)*dctmtx_8)+128;
    end
end

original_JPEG.JPEG_struct.coef_arrays=mat2cell(Quantized_DCT,height,width);
watermarked_JPEG.JPEG_struct.coef_arrays=mat2cell(watermarked_qc,height,width);

% stream_watermarked_DC = cell(1,length(watermarked_DC));
% for m2=1:length(watermarked_DC)
%     if m2==1
%         stream_watermarked_DC{m2}=huffman_dc(watermarked_DC(m2));
%     else
%         stream_watermarked_DC{m2}=huffman_dc(watermarked_DC(m2)-watermarked_DC(m2-1));
%     end
% end
% 
% 
% %% Huffman Coding
% %Original
% stream_DC=cell(1,length(zig_zag_coeff));
% stream_DC{1}=huffman_dc(zig_zag_coeff(1,1));
% stream_AC=cell(1,length(zig_zag_coeff));
% stream_AC{1}=huffman_ac(zig_zag_coeff(1,2:end));
% for m=2:length(zig_zag_coeff)
%     stream_DC{m}=huffman_dc(zig_zag_coeff(m,1)-zig_zag_coeff(m-1,1));
%     stream_AC{m}=huffman_ac(zig_zag_coeff(m,2:end));
% end

end

function [sorted_list,zig_zag_coeff]=prediction_list(TQ,DCTQ,order_zig,polarity)
p_reconstructed=zeros(size(TQ,1),size(TQ,2));
reconstructed=zeros(size(TQ,1),size(TQ,2));
zig_zag_coeff=zeros(size(TQ,1)/64*size(TQ,2),64);
DC_p=zeros(size(TQ,1)/64*size(TQ,2),1);
DC_p2=zeros(size(TQ,1)/64*size(TQ,2),1);
Dpcm=zeros(size(TQ,1)/64*size(TQ,2),1);
dctmtx_8=dctmtx(8);
sorted_list=zeros(size(TQ,1)/64*size(TQ,2),4);
P_TQ=zeros(size(TQ,1),size(TQ,2));
DC_prediction=zeros(size(TQ,1)/64*size(TQ,2),1);
for index_i=1:size(TQ,1)/8
    for index_j=1:size(TQ,2)/8
        i=1+8*(index_i-1);
        j=1+8*(index_j-1);
        reconstructed(i:i+7,j:j+7)=round(dctmtx_8'*(TQ(i:i+7,j:j+7).*DCTQ)*dctmtx_8)+128;
        k=index_j+(index_i-1)*size(TQ,2)/8;
        temp_tq=reshape(TQ(i:i+7,j:j+7)',1,64);
        zig_zag_coeff(k,:)=temp_tq(order_zig);
        P_TQ(i:i+7,j:j+7)=TQ(i:i+7,j:j+7);
        P_TQ(i,j)=0;
        p_reconstructed(i:i+7,j:j+7)=round(dctmtx_8'*(P_TQ(i:i+7,j:j+7).*DCTQ)*dctmtx_8)+128;
    end
end
k=0;
previous_DC=0;
for index_i=1:size(TQ,1)/8
    for index_j=1:size(TQ,2)/8
        i=1+8*(index_i-1);
        j=1+8*(index_j-1);
        current_qc=TQ(i:i+7,j:j+7);
        if mod(index_i+index_j,2)==polarity
            
            k=k+1;
            
            %First block
            if i==1 && j==1
                DC_prediction(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1);
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1);
                
                %Top right corner
            elseif i==1 && j==size(TQ,2)-7
                W=double(reconstructed(i:i+7,j-8:j-1));
                S=double(reconstructed(i+8:i+15,j:j+7));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_prediction(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); S(1,:)'-current_p(8,:)'])*8/DCTQ(1,1));
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                
                %Bottom left corner
            elseif i==size(TQ,1)-7 && j==1
                N=double(reconstructed(i-8:i-1,j:j+7));
                E=double(reconstructed(i:i+7,j+8:j+15));
                DC_prediction(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([N(8,:)'-current_p(1,:)'; E(:,1)-current_p(:,8)])*8/DCTQ(1,1));
                
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Bottom right corner
            elseif i==size(TQ,1)-7 && j==size(TQ,2)-7
                N=double(reconstructed(i-8:i-1,j:j+7));
                W=double(reconstructed(i:i+7,j-8:j-1));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_prediction(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); N(8,:)'-current_p(1,:)'])*8/DCTQ(1,1));
                
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Top line
            elseif i==1
                W=double(reconstructed(i:i+7,j-8:j-1));
                E=double(reconstructed(i:i+7,j+8:j+15));
                S=double(reconstructed(i+8:i+15,j:j+7));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_prediction(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); E(:,1)-current_p(:,8); S(1,:)'-current_p(8,:)'])*8/DCTQ(1,1));
                
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Left line
            elseif j==1
                N=double(reconstructed(i-8:i-1,j:j+7));
                E=double(reconstructed(i:i+7,j+8:j+15));
                S=double(reconstructed(i+8:i+15,j:j+7));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_prediction(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([N(8,:)'-current_p(1,:)'; E(:,1)-current_p(:,8); S(1,:)'-current_p(8,:)'])*8/DCTQ(1,1));
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Right line
            elseif j== size(TQ,2)-7
                N=double(reconstructed(i-8:i-1,j:j+7));
                W=double(reconstructed(i:i+7,j-8:j-1));
                current_p=p_reconstructed(i:i+7,j:j+7);
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                DC_prediction(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); N(8,:)'-current_p(1,:)'])*8/DCTQ(1,1));
                
                %Bottom line
            elseif i== size(TQ,1)-7
                N=double(reconstructed(i-8:i-1,j:j+7));
                W=double(reconstructed(i:i+7,j-8:j-1));
                E=double(reconstructed(i:i+7,j+8:j+15));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_prediction(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); N(8,:)'-current_p(1,:)'; E(:,1)-current_p(:,8)])*8/DCTQ(1,1));
                
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Rest
            else
                N=double(reconstructed(i-8:i-1,j:j+7));
                W=double(reconstructed(i:i+7,j-8:j-1));
                E=double(reconstructed(i:i+7,j+8:j+15));
                S=double(reconstructed(i+8:i+15,j:j+7));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_prediction(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); N(8,:)'-current_p(1,:)'; E(:,1)-current_p(:,8); S(1,:)'-current_p(8,:)'])*8/DCTQ(1,1));
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
            end
            
            sorted_list(k,:)=[i,j,DC_prediction(index_j+(index_i-1)*size(TQ,2)/8),sum(zig_zag_coeff(k,2:end)==0)];
        end
        previous_DC=current_qc(1,1);
    end
end
sorted_list=sorted_list(1:k,:);
sorted_list=sortrows(sorted_list,-4);
end

