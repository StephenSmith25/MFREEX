function [X,Y,C_2D] = temperature_profile(temperature, time)

    nodeN = 81;
    folder = ['./Preform Heating/Oil Bath/Heating_cycle_',num2str(temperature),'_tc5','/'];
    for i = 1:1:151
        T_profile(i,1) = i;
        for j = 1:1:5
            name = ['NT11_node_',num2str((i-1)*5+j),'.csv'];
            temp = csvread([folder,name]);
            T_profile(i,j+1) = interp1(temp(:,1),temp(:,2),time);
        end
        T_profile(i,7) = sum(T_profile(i,2:6))/5;
    end
    inter_x = [1:150/(nodeN-1):151];
    inter_x = inter_x';
    for j = 1:1:6
        interp_T_profile(:,j) = interp1(T_profile(:,1),T_profile(:,j+1),inter_x);
    end
    result = interp_T_profile;
    %   draw figures
    load('COORDs.mat')
    
    for i = 1:1:151
        for j = 1:1:5
            X(i,j)=raw((i-1)*5+j,2);
            Y(i,j)=raw((i-1)*5+j,3);
            C_1D(i,j) = T_profile(i,7);
            C_2D(i,j) = T_profile(i,j+1);
            C_0D = mean(T_profile);
        end
    end
% 
%     h1 = figure(1);
%     h1 = pcolor(X,Y,C_1D);
%     set(h1, 'LineStyle','none')
%     %axis off;
%     axis([0,80,-80,0]);
%     %%caxis([101 102])
%     colorbar;
% 
%     h2 = figure(2);
%     h2 = pcolor(X,Y,C_2D);
%     set(h2, 'LineStyle','none')
%     %axis off;
%     axis([0,80,-80,0]);
%     %%caxis([101 105])
%     colorbar;
  
end
