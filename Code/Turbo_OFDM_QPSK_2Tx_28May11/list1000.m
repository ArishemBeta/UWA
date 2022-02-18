%Created by Jun Tao
%
function [packet_name, Sync]=list1000()

packet_name(1,:)= '3011553F010_C0_S6'; Sync(1,:)=[5,0];
packet_name(2,:)= '3011554F010_C0_S5'; Sync(2,:)=[-5,0];

packet_name(3,:)= '3011753F010_C0_S6'; Sync(3,:)=[5,0];
packet_name(4,:)= '3011754F010_C0_S5'; Sync(4,:)=[-25,0]; %ofdm, worse

packet_name(5,:)= '3011953F010_C0_S6'; Sync(5,:)=[5,0];
packet_name(6,:)= '3011954F010_C0_S5'; Sync(6,:)=[-5,0];

packet_name(7,:)= '3012153F010_C0_S6'; Sync(7,:)=[5,0];
packet_name(8,:)= '3012154F010_C0_S5'; Sync(8,:)=[-8,0];

packet_name(9,:)= '3012353F010_C0_S6'; Sync(9,:)=[5,0];
packet_name(10,:)= '3012354F010_C0_S5'; Sync(10,:)=[-5,0];

packet_name(11,:)= '3020156F010_C0_S6'; Sync(11,:)=[6,0];

packet_name(12,:)= '3020353F010_C0_S6'; Sync(12,:)=[0,0];

packet_name(13,:)= '3020553F010_C0_S6'; Sync(13,:)=[5,0];

packet_name(14,:)= '3020753F010_C0_S6'; Sync(14,:)=[0,0];

packet_name(15,:)= '3020953F010_C0_S6'; Sync(15,:)=[5,0];

%------end of data loading------
