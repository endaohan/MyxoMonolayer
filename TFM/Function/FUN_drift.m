function D_Glb = FUN_drift(route,routeOut,Parameter)

% For debugging
%{
% nCh = 2; 
routeDrift = [folder '\' fileName];
route = routeDrift; 
routeOut = [route '\Drift']; 
%}

nCh = Parameter.nCh-1; 

% Low resolution displacement
% col 1: raw drift data; col 2: smoothed drift data
data = cell(nCh,2); 
% Mean data
route_LowRes = [route '\Beads_1\LowRes_Data']; 
PIVList_LR = FUN_FileName(route_LowRes,'.mat','PIV_','No'); 
N_LR = length(PIVList_LR); 
u_mean = zeros(N_LR,nCh); 
v_mean = zeros(N_LR,nCh); 

for ind = 1:1:nCh
    % Set input routes
    route_LowRes = [route '\Beads_' num2str(ind) '\LowRes_Data']; 

    %{
    route_LowRes = [folder '\' fileList{ind_file} '\Beads_1\PIV_LowRes_Data']; 
    route_TFM = [folder '\' fileList{ind_file} '\Beads_1\PIV_Data']; 
    %}

    % List of file names
    PIVList_LR = FUN_FileName(route_LowRes,'.mat','PIV_','No'); 
    N_LR = length(PIVList_LR); 


    %% Get information from the first frame
    % Low resolution data for drift
    V_LR = zeros(N_LR,2); 
    % T = linspace(1, N_LR, N_LR)*Parameter.tscl/60;    % Time, unit: min
    T = linspace(1, N_LR, N_LR);    % Time, unit: frame 


    %% Drift
    for i = 1:1:N_LR
        load([route_LowRes '\' PIVList_LR{i} '.mat'],'Vectors'); 

        u = rmoutliers(Vectors.u); 
        v = rmoutliers(Vectors.v);  

        % Four colomns: mean u, mean v, stdev u, stdev v. 
        V_LR(i,1) = mean(u); 
        V_LR(i,2) = mean(v); 
        V_LR(i,3) = std(u); 
        V_LR(i,4) = std(v); 
    end
    clear Vectors
    
    % Remove the zero frequency mode
    V_LR(:,1) = V_LR(:,1)-mean(V_LR(:,1)); 
    V_LR(:,2) = V_LR(:,2)-mean(V_LR(:,2)); 
    
    % Smooth the drift
    %{
    if N_LR >= 10
        V_smooth(:,1) = FUN_GaussSmooth(V_LR(:,1),10); 
        V_smooth(:,2) = FUN_GaussSmooth(V_LR(:,2),10); 
    else
        V_smooth(:,1) = FUN_GaussSmooth(V_LR(:,1),5); 
        V_smooth(:,2) = FUN_GaussSmooth(V_LR(:,2),5);
    end
    %}
    
    % Mean data matrix
    u_mean(:,ind) = V_LR(:,1); 
    v_mean(:,ind) = V_LR(:,2); 
    
    % Put the data in the cell
    data{ind,1} = V_LR; 
    % data{ind,2} = V_smooth; 
    
end

% mean displacement
U = mean(u_mean,2); 
V = mean(v_mean,2); 


% Make plot
hd_drift = figure('Unit','Inch','Position',[2 2 6 3]); 
subplot(1,2,1)
hold on
for j = 1:1:nCh
    errorbar(T,data{j,1}(:,1),data{j,1}(:,3),'b-')
    % plot(T,data{j,2}(:,1)*1E6,'b--')
end
plot(T,U,'r-')
hold off
xlabel('Time (frame)')
ylabel('Drift in x (pixel)')
% legend('Drift in x','Drift in y','Location','NorthEast')
% axis([0 60 -0.2 1])
box on

subplot(1,2,2)
hold on
for j = 1:1:nCh
    errorbar(T,data{j,1}(:,2),data{j,1}(:,4),'b-')
    % plot(T,data{j,2}(:,2)*1E6,'b--')
end
plot(T,V,'r-')
hold off
xlabel('Time (frame)')
ylabel('Drift in y (pixel)')
box on

% Save the results
% D_Glb.V_LR = V_LR; 
D_Glb.V_LR = data; 
D_Glb.V_mean = [U V]; 
% D_Glb.V_sm = V_smooth; 
saveas(hd_drift, [routeOut '\GlobalDrift.png']);
save([routeOut '\GlobalDrift.mat'],'D_Glb'); 


end


