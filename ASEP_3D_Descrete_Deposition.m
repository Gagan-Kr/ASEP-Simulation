%% Lattice Config

LattS = [ 10 250 100 ];   %   [ Height Length Width ]

Latt = false(LattS);
J = false(cat(2,LattS,6));
FlowTimeCount = 0;

f1 = figure(1);
[x,y,z] = ind2sub(size(Latt),find(Latt));
PlotLatt = scatter3(x,y,z,2,"filled");
xlabel("Height"); ylabel("Length"); zlabel("Weidth");
PlottLattTitle = title("Gas flow","(Time "+num2str(FlowTimeCount)+")");
xlim([0 LattS(1)]); ylim([0 LattS(2)]); zlim([0 LattS(3)]);
view(100.7794,37.6074);

figure(3)
SidViewFlowH = heatmap(sum(Latt,3),"Colormap",jet);
SidViewFlowH.Title = ["Side View","(Time "+ num2str(FlowTimeCount)+")"];
grid off;

figure(4)
TopViewFlowH = heatmap(rot90(squeeze(sum(Latt,1)),1),"Colormap",jet);
TopViewFlowH.Title = ["Top View","(Time "+ num2str(FlowTimeCount)+")"];
grid off;

%% Substrate

[SubsDim2, SubsDim3] =  deal( 30:LattS(2) , 1:LattS(3) );

Substrate = zeros(1,LattS(2),LattS(3));
Substrate(1,LattS(2),LattS(3)) = 0.001;
SubsTimeCount = 0;

f2 = figure(2);
SubstrateH = heatmap(rot90(squeeze(Substrate),1),"Colormap",turbo);
grid off;

figure(5)
SubstrateHPlot = plot(1:LattS(2),zeros(1,LattS(2)));
grid on;

figure(6)
SubstrateLPlot = plot(1:LattS(3),zeros(1,LattS(3)));
grid on;

%% Injection & Ejection window

[IL3, IL1, EL3, EL1] = deal( [] , [] , [] , [] );
[IR3, IR1, ER3, ER1] = deal( [] , [] , 1:LattS(3) , 1:LattS(1) );
[IU3, IU2, EU3, EU2] = deal( [] , [] , [] , [] );
[ID3, ID2, ED3, ED2] = deal( 36:65 , 11:40 , [] , [] );

%% Simukation Settngs

Alfa = 0.03;
Beta = 0.1;

[ XBias , YBias , ZBias ]  = deal( 0.07 , 0 , 0 );

Deposition = true; RxnRate = 1;

Simulation = true; Speed = 100;
Time = 5000;

[p, q] = deal( 1/6+(XBias/6) , 1/6-(XBias/6) );
[r, s] = deal( 1/6+(YBias/6) , 1/6-(YBias/6) );
[u, v] = deal( 1/6+(ZBias/6) , 1/6-(ZBias/6) );

%% Run Simulation

figure(1);
tic
for t = 1:Time
    FlowTimeCount = FlowTimeCount + 1;

    % Injection & Ejection

    Latt(IR1,end,IR3) = Latt(IR1,end,IR3) | rand(length(IR1),1,length(IR3)) < Alfa;
    Latt(end,ID2,ID3) = Latt(end,ID2,ID3) | rand(1,length(ID2),length(ID3)) < Alfa;
    Latt(IL1,1,IL3) = Latt(IL1,1,IL3) | rand(length(IL1),1,length(IL3)) < Alfa;
    Latt(1,IU2,IU3) = Latt(1,IU2,IU3) | rand(1,length(IU2),length(IU3)) < Alfa;

    Latt(ER1,end,ER3) = Latt(ER1,end,ER3) & rand(length(ER1),1,length(ER3)) > Beta;
    Latt(end,ED2,ED3) = Latt(end,ED2,ED3) & rand(1,length(ED2),length(ED3)) > Beta;
    Latt(EL1,1,EL3) = Latt(EL1,1,EL3) & rand(length(EL1),1,length(EL3)) > Beta;
    Latt(1,EU2,EU3) = Latt(1,EU2,EU3) & rand(1,length(EU2),length(EU3)) > Beta;

    % Asymetric Exclusion Process

    R = rand(LattS);
    R(Latt == false) = nan;
    M = randperm(6);

    % Deposition

    if Deposition && rand < RxnRate
        SubsTimeCount = SubsTimeCount + 1;
        JD = R(1,SubsDim2,SubsDim3)>=p+q & R(1,SubsDim2,SubsDim3)<p+q+r;
        Substrate(1,SubsDim2,SubsDim3) = Substrate(1,SubsDim2,SubsDim3) + JD;
    end

    for i = 1:6

        if M(i) == 1   % Right hop
            J(:,1:end-1,:,1) = R(:,1:end-1,:)<p & ~Latt(:,2:end,:);
            Latt(circshift(J(:,:,:,1),1,2)) = true;
        end

        if M(i) == 2   % Left hop
            J(:,2:end,:,2) = R(:,2:end,:)>=p & R(:,2:end,:)<p+q & ~Latt(:,1:end-1,:);
            Latt(circshift(J(:,:,:,2),-1,2)) = true;
        end
 
        if M(i) == 3   % Up hop
            J(:,:,1:end-1,3) = (R(:,:,1:end-1)>=p+q & ~Latt(:,:,2:end) & R(:,:,1:end-1)<p+q+r);
            Latt(circshift(J(:,:,:,3),1,3)) = true;
        end

        if M(i) == 4   % Down hop
            J(:,:,2:end,4) = R(:,:,2:end)>=p+q+r & ~Latt(:,:,1:end-1) & R(:,:,2:end)<p+q+r+s;
            Latt(circshift(J(:,:,:,4),-1,3)) = true;
        end

        if M(i) == 5   % Inward hop
            J(1:end-1,:,:,5) = R(1:end-1,:,:)>=p+q+r+s & ~Latt(2:end,:,:) & R(1:end-1,:,:)<p+q+r+s+u;
            Latt(circshift(J(:,:,:,5),1,1)) = true;
        end

        if M(i) == 6   % Outward hop
            J(2:end,:,:,6) = R(2:end,:,:)>=p+q+r+s+u & ~Latt(1:end-1,:,:);
            Latt(circshift(J(:,:,:,6),-1,1)) = true;
        end
    end

    Latt(J(:,:,:,1)|J(:,:,:,2)|J(:,:,:,3)|J(:,:,:,4)|J(:,:,:,5)|J(:,:,:,6)) = false;

    % Simulaton

    if Simulation && floor(FlowTimeCount/Speed) == FlowTimeCount/Speed
        [x,y,z] = ind2sub(size(Latt),find(Latt));
        set(PlotLatt,'XData',x,'YData',y,'ZData',z)
        title("Gas flow","(Time "+num2str(FlowTimeCount)+")");
        SubstrateData = rot90(squeeze(Substrate));
        set(SubstrateH,'ColorData',SubstrateData);
        SubstrateH.Title = ["Deposition","(Time "+num2str(SubsTimeCount)+")"];
        SidViewFlowH.ColorData = sum(Latt,3);
        TopViewFlowH.ColorData = rot90(squeeze(sum(Latt,1)),1);
        SidViewFlowH.Title = ["Side View","(Time "+ num2str(FlowTimeCount)+")"];
        TopViewFlowH.Title = [ "Top View","(Time "+ num2str(FlowTimeCount)+")"];
        Substrate2dH = sum(SubstrateData((LattS(3)/2)-1:LattS(3)/2+2,:),1);
        [ ~ , Ind] = max(Substrate2dH);
        if Ind < 2
           Ind = 2;
        end
        Substrate2dV = sum(SubstrateData(:,Ind-1:Ind+2),2);
        set(SubstrateHPlot,'YData',Substrate2dH/4);
        set(SubstrateLPlot,'YData',Substrate2dV/4);

        drawnow;
    end

    if floor((t/Time)*100) == (t/Time)*100
        clc; fprintf('Simulating... (%d%%)\n', (t/Time)*100);
    end
end

toc
