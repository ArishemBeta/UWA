%Created by Jun Tao
%

function [packet_name, Sync]=list200()

packet_name(1,:)= '3011554F010_C0_S3'; Sync(1,:)=[5,0]; %-ofdm
packet_name(2,:)= '3011554F010_C0_S4'; Sync(2,:)=[5,0]; %ofdm
packet_name(3,:)= '3011554F010_C1_S3'; Sync(3,:)=[5,0]; %-ofdm (bad)

packet_name(4,:)= '3011754F010_C0_S3'; Sync(4,:)=[5,0]; %-ofdm
packet_name(5,:)= '3011754F010_C0_S4'; Sync(5,:)=[5,0]; %ofdm
packet_name(6,:)= '3011754F010_C1_S3'; Sync(6,:)=[5,0]; %-ofdm

packet_name(7,:)= '3011954F010_C0_S3'; Sync(7,:)=[10,0];
packet_name(8,:)= '3011953F010_C0_S4'; Sync(8,:)=[0,0]; %ofdm
packet_name(9,:)= '3011954F010_C1_S3'; Sync(9,:)=[0,0];

packet_name(10,:)= '3012154F010_C0_S3'; Sync(10,:)=[10,0];
packet_name(11,:)= '3012154F010_C0_S4'; Sync(11,:)=[7,0]; %ofdm
packet_name(12,:)= '3012154F010_C1_S3'; Sync(12,:)=[-8,0];

packet_name(13,:)= '3012354F010_C0_S3'; Sync(13,:)=[10,0];
packet_name(14,:)= '3012353F010_C0_S4'; Sync(14,:)=[7,0]; %ofdm
packet_name(15,:)= '3012354F010_C1_S3'; Sync(15,:)=[10,0];

packet_name(16,:)= '3020157F010_C0_S3'; Sync(16,:)=[10,0];
packet_name(17,:)= '3020157F010_C0_S4'; Sync(17,:)=[7,0]; %ofdm
packet_name(18,:)= '3020157F010_C1_S3'; Sync(18,:)=[10,0];

packet_name(19,:)= '3020354F010_C0_S3'; Sync(19,:)=[10,0];
packet_name(20,:)= '3020354F010_C0_S4'; Sync(20,:)=[7,0]; %ofdm
packet_name(21,:)= '3020354F010_C1_S3'; Sync(21,:)=[10,0];

packet_name(22,:)= '3020554F010_C0_S3'; Sync(22,:)=[10,0];
packet_name(23,:)= '3020554F010_C0_S4'; Sync(23,:)=[-10,0]; %ofdm
packet_name(24,:)= '3020554F010_C1_S3'; Sync(24,:)=[10,0];

packet_name(25,:)= '3020754F010_C0_S3'; Sync(25,:)=[10,0];
packet_name(26,:)= '3020754F010_C0_S4'; Sync(26,:)=[10,0]; %-ofdm
packet_name(27,:)= '3020754F010_C1_S3'; Sync(27,:)=[12,0];

packet_name(28,:)= '3020954F010_C0_S3'; Sync(28,:)=[12,0];
packet_name(29,:)= '3020954F010_C0_S4'; Sync(29,:)=[10,0]; %-ofdm
packet_name(30,:)= '3020954F010_C1_S3'; Sync(30,:)=[12,0];

packet_name(31,:)= '3021154F010_C0_S3'; Sync(31,:)=[12,0];
packet_name(32,:)= '3021153F010_C0_S4'; Sync(32,:)=[10,0]; %-ofdm
packet_name(33,:)= '3021154F010_C1_S3'; Sync(33,:)=[-10,0];

packet_name(34,:)= '3021354F010_C1_S3'; Sync(34,:)=[12,0];
packet_name(35,:)= '3021353F010_C0_S4'; Sync(35,:)=[5,0];
packet_name(36,:)= '3021354F010_C0_S3'; Sync(36,:)=[12,0];

packet_name(37,:)= '3021554F010_C0_S3'; Sync(37,:)=[10,0];
packet_name(38,:)= '3021554F010_C0_S4'; Sync(38,:)=[5,0];
packet_name(39,:)= '3021554F010_C1_S3'; Sync(39,:)=[10,0];

packet_name(40,:)= '3021754F010_C0_S3'; Sync(40,:)=[12,0];
packet_name(41,:)= '3021754F010_C0_S4'; Sync(41,:)=[5,0];
packet_name(42,:)= '3021754F010_C1_S3'; Sync(42,:)=[12,0];

packet_name(43,:)= '3021954F010_C0_S3'; Sync(43,:)=[10,0];
packet_name(44,:)= '3021953F010_C0_S4'; Sync(44,:)=[5,0];
packet_name(45,:)= '3021954F010_C1_S3'; Sync(45,:)=[10,0];

