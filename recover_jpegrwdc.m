function [recovered_payload, recovered_qc] = recover_jpegrwdc(watermarked_qc,q_factor)
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

[col row]=size(image);
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

[sorted_list_odd,~]=prediction_list(watermarked_qc,DCTQ,order_zig,1);
%recover
[col, row]=size(watermarked_qc);

flag=0;
payload_length_limit=ceil(log2(col*row/64));
recovered_side_information_payload_size=ones(1,payload_length_limit);
recovered_payload_limit=bi2de(recovered_side_information_payload_size);
side_information_length=0;
recovered_qc=watermarked_qc;
payload_size=0;
recovered_payload=zeros(1,2^payload_length_limit);

for i=1:length(sorted_list_odd)
    if flag==0
        dc_prediction=sorted_list_odd(i,3);
        i1=sorted_list_odd(i,1);
        j1=sorted_list_odd(i,2);
        if dc_prediction == 0 || dc_prediction == -1
            if side_information_length < payload_length_limit
                side_information_length=side_information_length+1;
                recovered_side_information_payload_size(side_information_length)=0;
                if side_information_length ==payload_length_limit
                    recovered_payload_limit=bi2de(recovered_side_information_payload_size);
                end
            else
                payload_size=payload_size+1;
                recovered_payload(payload_size)=0;
                if payload_size == recovered_payload_limit/2
                    flag=1;
                end
            end
        elseif dc_prediction == 1 || dc_prediction == -2
            if side_information_length < payload_length_limit
                side_information_length=side_information_length+1;
                recovered_side_information_payload_size(side_information_length)=1;
                if side_information_length ==payload_length_limit
                    recovered_payload_limit=bi2de(recovered_side_information_payload_size);
                end
            else
                payload_size=payload_size+1;
                recovered_payload(payload_size)=1;
                if payload_size == recovered_payload_limit/2
                    flag=1;
                end
            end
            if dc_prediction == 1
                recovered_qc(i1,j1)=watermarked_qc(i1,j1)-1;
            elseif dc_prediction == -2
                recovered_qc(i1,j1)=watermarked_qc(i1,j1)+1;
            end
        elseif dc_prediction < -1
            recovered_qc(i1,j1)=watermarked_qc(i1,j1)+1;
        elseif dc_prediction > 1
            recovered_qc(i1,j1)=watermarked_qc(i1,j1)-1;
        end
    end
    if  i == length(sorted_list_odd) && flag == 0
        flag = 2;
        pause
    end
end
recovered_payload(payload_size+1:end)=[];
%polarity == 0
[sorted_list,~]=prediction_list(recovered_qc,DCTQ,order_zig,0);

flag=0;
%embed even set
for i=1:length(sorted_list)
    if flag==0
        dc_prediction=sorted_list(i,3);
        i1=sorted_list(i,1);
        j1=sorted_list(i,2);
        if dc_prediction == 0 || dc_prediction == -1
            payload_size=payload_size+1;
            recovered_payload(payload_size)=0;
            if payload_size == recovered_payload_limit
                flag=1;
            end
        elseif dc_prediction == 1 || dc_prediction == -2
            payload_size=payload_size+1;
            recovered_payload(payload_size)=1;
            if payload_size == recovered_payload_limit
                flag=1;
            end
            if dc_prediction == 1
                recovered_qc(i1,j1)=watermarked_qc(i1,j1)-1;
            elseif dc_prediction == -2
                recovered_qc(i1,j1)=watermarked_qc(i1,j1)+1;
            end
        elseif dc_prediction < -1
            recovered_qc(i1,j1)=watermarked_qc(i1,j1)+1;
        elseif dc_prediction > 1
            recovered_qc(i1,j1)=watermarked_qc(i1,j1)-1;
        end
    end
end


temp_payload=[recovered_payload(recovered_payload_limit/2+1:end) recovered_payload(1:recovered_payload_limit/2)];
recovered_payload=temp_payload;

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
DC_p3=zeros(size(TQ,1)/64*size(TQ,2),1);
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
                DC_p3(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1);
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1);
                
                %Top right corner
            elseif i==1 && j==size(TQ,2)-7
                W=double(reconstructed(i:i+7,j-8:j-1));
                S=double(reconstructed(i+8:i+15,j:j+7));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_p3(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); S(1,:)'-current_p(8,:)'])*8/DCTQ(1,1));
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                
                %Bottom left corner
            elseif i==size(TQ,1)-7 && j==1
                N=double(reconstructed(i-8:i-1,j:j+7));
                E=double(reconstructed(i:i+7,j+8:j+15));
                DC_p3(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([N(8,:)'-current_p(1,:)'; E(:,1)-current_p(:,8)])*8/DCTQ(1,1));
                
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Bottom right corner
            elseif i==size(TQ,1)-7 && j==size(TQ,2)-7
                N=double(reconstructed(i-8:i-1,j:j+7));
                W=double(reconstructed(i:i+7,j-8:j-1));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_p3(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); N(8,:)'-current_p(1,:)'])*8/DCTQ(1,1));
                
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Top line
            elseif i==1
                W=double(reconstructed(i:i+7,j-8:j-1));
                E=double(reconstructed(i:i+7,j+8:j+15));
                S=double(reconstructed(i+8:i+15,j:j+7));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_p3(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); E(:,1)-current_p(:,8); S(1,:)'-current_p(8,:)'])*8/DCTQ(1,1));
                
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Left line
            elseif j==1
                N=double(reconstructed(i-8:i-1,j:j+7));
                E=double(reconstructed(i:i+7,j+8:j+15));
                S=double(reconstructed(i+8:i+15,j:j+7));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_p3(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([N(8,:)'-current_p(1,:)'; E(:,1)-current_p(:,8); S(1,:)'-current_p(8,:)'])*8/DCTQ(1,1));
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Right line
            elseif j== size(TQ,2)-7
                N=double(reconstructed(i-8:i-1,j:j+7));
                W=double(reconstructed(i:i+7,j-8:j-1));
                current_p=p_reconstructed(i:i+7,j:j+7);
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                DC_p3(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); N(8,:)'-current_p(1,:)'])*8/DCTQ(1,1));
                
                %Bottom line
            elseif i== size(TQ,1)-7
                N=double(reconstructed(i-8:i-1,j:j+7));
                W=double(reconstructed(i:i+7,j-8:j-1));
                E=double(reconstructed(i:i+7,j+8:j+15));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_p3(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); N(8,:)'-current_p(1,:)'; E(:,1)-current_p(:,8)])*8/DCTQ(1,1));
                
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
                %Rest
            else
                N=double(reconstructed(i-8:i-1,j:j+7));
                W=double(reconstructed(i:i+7,j-8:j-1));
                E=double(reconstructed(i:i+7,j+8:j+15));
                S=double(reconstructed(i+8:i+15,j:j+7));
                current_p=p_reconstructed(i:i+7,j:j+7);
                DC_p3(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-round(mean([W(:,8)-current_p(:,1); N(8,:)'-current_p(1,:)'; E(:,1)-current_p(:,8); S(1,:)'-current_p(8,:)'])*8/DCTQ(1,1));
                Dpcm(index_j+(index_i-1)*size(TQ,2)/8)=current_qc(1,1)-previous_DC;
            end
            
            sorted_list(k,:)=[i,j,DC_p3(index_j+(index_i-1)*size(TQ,2)/8),sum(zig_zag_coeff(k,2:end)==0)];
        end
        previous_DC=current_qc(1,1);
    end
end
sorted_list=sorted_list(1:k,:);
sorted_list=sortrows(sorted_list,-4);
% sorted_list2=sortrows(sorted_list2,-4);
% sorted_list3=sortrows(sorted_list3,-4);

end
