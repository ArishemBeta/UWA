clear all; close all

N_chnn=4;

for k=0: -5: 0
    k
    step1=0;
    %Find_Dpr_Psbd2Bsbd_0822(step1,k);
    %Find_Dpr_Psbd2Bsbd_0837(step1,k);
    %Find_Dpr_Psbd2Bsbd_1007(step1,k);
    %Find_Dpr_Psbd2Bsbd_1123(step1,k, N_chnn);
    %Find_Dpr_Psbd2Bsbd_0708(step1,k, N_chnn);
    Find_Dpr_Psbd2Bsbd_0853(step1,k, N_chnn);
    pause(2); %delay 2s
    main
%     fr2=fopen('sync.dat','a');
%     fprintf(fr2,'%d %d %f\n', k, Error_Num, BER);
%     fclose(fr2); 
end