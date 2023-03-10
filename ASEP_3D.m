%% Lattice Config

%         Height Length Width
LattSize = [ 50 , 100 , 70 ];

Latt = zeros(LattSize+2);
Latt(1,:,:) = 1; Latt(end,:,:) = 1;
Latt(:,1,:) = 1; Latt(:,end,:) = 1;
Latt(:,:,1) = 1; Latt(:,:,end) = 1;

FlowTimeCount = 0;
 
%% Injection & Ejection window

[InjL3, InjL1, EjeL3, EjeL1] = deal( 16:35 , 16:35 , [] , [] );
[InjR3, InjR1, EjeR3, EjeR1] = deal( [] , [] , 16:35 , 16:35 );
[InjU3, InjU2, EjeU3, EjeU2] = deal( [] , [] , [] , [] );
[InjD3, InjD2, EjeD3, EjeD2] = deal( [] , [] , [] , [] );

AlfaR = ~isempty(InjR1) && ~isempty(InjR3);
AlfaL = ~isempty(InjL1) && ~isempty(InjL3);
AlfaU = ~isempty(InjU2) && ~isempty(InjU3);
AlfaD = ~isempty(InjD3) && ~isempty(InjD2);

BetaR = ~isempty(EjeR1) && ~isempty(EjeR3);
BetaL = ~isempty(EjeL1) && ~isempty(EjeL3);
BetaU = ~isempty(EjeU3) && ~isempty(EjeU2);
BetaD = ~isempty(EjeD3) && ~isempty(EjeD2);

% LattB(InjR1+1,end-1,InjR3+1) = 1;
% LattB(end-1,InjD2+1,InjD3+1) = 1;
% LattB(InjL1+1,2,InjL3+1) = 1;
% LattB(2,InjU2+1,InjU3+1) = 1;
% 
% LattB(EjeR1+1,end-1,EjeR3+1) = -1;
% LattB(end-1,EjeD2+1,EjeD3+1) = -1;
% LattB(EjeL1+1,2,EjeL3+1) = -1;
% LattB(2,EjeU2+1,EjeL3+1) = -1;

figure(1)
% TopViewFlowH = heatmap(rot90(squeeze(sum(Latt(2:end-1,2:end-1,2:end-1),1)),1),"Colormap",jet);
% TopViewFlowH = spy(rot90(squeeze(sum(Latt(2:end-1,2:end-1,2:end-1),1)),1));
% TopViewFlowH.Title = ["Top View","(Time "+ num2str(FlowTimeCount)+")"];
%grid off;
[x,y,z]=ind2sub(size(Latt(2:end-1,2:end-1,2:end-1)),find(Latt(2:end-1,2:end-1,2:end-1)));
PlotLatt = scatter3(x,y,z);
xlim([0 LattSize(1)]);
ylim([0 LattSize(2)]);
zlim([0 LattSize(3)]);

figure(2)
SidViewFlowH = heatmap(sum(Latt(2:end-1,2:end-1,2:end-1),3),"Colormap",jet);
SidViewFlowH.Title = ["Side View","(Time "+ num2str(FlowTimeCount)+")"];
grid off;



%% Substrate

[SubsDim2, SubsDim3, SubstrtHeight] =  deal( 1:LattSize(2) , 1:LattSize(3) , 0 );

Substrate = zeros(1,LattSize(2),LattSize(3));
Substrate(1,LattSize(2),LattSize(3)) = 0.001;

figure(3)
SubstrateH = heatmap(rot90(squeeze(Substrate),1),"Colormap",turbo);
grid off;

SubsTimeCount = 0;

%% Simukation Settngs

Alfa = 1;
Beta = 1;

XBias = 0;
YBias = 0;
ZBias = 0;

[p, q] = deal( 1/6+(XBias/6) , 1/6-(XBias/6) );
[r, s] = deal( 1/6+(YBias/6) , 1/6-(YBias/6) );
[u, v] = deal( 1/6+(ZBias/6) , 1/6-(ZBias/6) );

Deposition = false;

Simulation = true; Speed = 1;

Time = 1000;


%% Run Simulation
tic 
for t = 1:Time

    FlowTimeCount = FlowTimeCount + 1;

    % Injection

    if AlfaR
        AlfaRateR = Alfa*(1-Latt(InjR1+1,end-1,InjR3+1));
        RandInjR = rand(size(AlfaRateR));
        Latt(InjR1+1,end-1,InjR3+1) = Latt(InjR1+1,end-1,InjR3+1) + (RandInjR<AlfaRateR);
    end

    if AlfaD
        AlfaRateD =  Alfa*(1-Latt(end-1,InjD2+1,InjD3+1));
        RandInjD = rand(size(AlfaRateD));
        Latt(end-1,InjD2+1,InjD3+1) = Latt(end-1,InjD2+1,InjD3+1) + (RandInjD<AlfaRateD);
    end

    if AlfaL
        AlfaRateL = Alfa*(1-Latt(InjL1+1,2,InjL3+1));
        RandInjL = rand(size(AlfaRateL));
        Latt(InjL1+1,2,InjL3+1) = Latt(InjL1+1,2,InjL3+1) + (RandInjL<AlfaRateL);
    end

    if AlfaU
        AlfaRateU =  Alfa*(1-Latt(2,InjU2+1,InjU3+1));
        RandInjU = rand(size(AlfaRateU));
        Latt(2,InjU2+1,InjU3+1) = Latt(2,InjU2+1,InjU3+1) + (RandInjU<AlfaRateU);
    end


    % Ejection

    if BetaR
        BetaRateR =  Beta*Latt(EjeR1+1,end-1,EjeR3+1);
        RandEjeR = rand(size(BetaRateR));
        Latt(EjeR1+1,end-1,EjeR3+1) = Latt(EjeR1+1,end-1,EjeR3+1) - (RandEjeR<BetaRateR);
    end

    if BetaD
        BetaRateD =  Beta*Latt(end-1,EjeD2+1,EjeD3+1);
        RandEjeD = rand(size(BetaRateD));
        Latt(end-1,EjeD2+1,EjeD3+1) = Latt(end-1,EjeD2+1,EjeD3+1) - (RandEjeD<BetaRateD);
    end

    if BetaL
        BetaRateL =  Beta*Latt(EjeL1+1,2,EjeL3+1);
        RandEjeL = rand(size(BetaRateL));
        Latt(EjeL1+1,2,EjeL3+1) = Latt(EjeL1+1,2,EjeL3+1) - (RandEjeL<BetaRateL);
    end

    if BetaU
        BetaRateU =  Betaa*Latt(2,EjeU2+1,EjeL3+1);
        RandEjeU = rand(size(BetaRateU));
        Latt(2,EjeU2+1,EjeU3+1) = Latt(2,EjeU2+1,EjeU3+1) - (RandEjeU<BetaRateU);
    end

    % Asymetric Exclusion Process

    Rand = rand(LattSize);

    M = randperm(6);

    for i = 1:6

        if M(i) == 1
            if i ~= 1
                J = Rand >= R;
                R = R + p*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,3:end,2:end-1)));
                J = J & Rand < R;
            else
                R = p*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,3:end,2:end-1)));
                J = Rand < R;
            end
            Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - J;
            Latt(2:end-1,3:end-1,2:end-1) = Latt(2:end-1,3:end-1,2:end-1) + J(:,1:end-1,:);
        end

        if M(i) == 2
            if i ~= 1
                J = Rand >= R;
                R = R + q*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,1:end-2,2:end-1)));
                J = J & Rand < R;
            else
                R = q*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,1:end-2,2:end-1)));
                J = Rand < R;
            end
            Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - J;
            Latt(2:end-1,2:end-2,2:end-1) = Latt(2:end-1,2:end-2,2:end-1) + J(:,2:end,:);
        end
 
        if M(i) == 3
            if i ~= 1
                J = Rand >= R;
                R = R + r*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,2:end-1,3:end)));
                J = J & Rand < R;
            else
                R = r*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,2:end-1,3:end)));
                J = Rand < R;
            end
            Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - J;
            Latt(2:end-1,2:end-1,3:end-1) = Latt(2:end-1,2:end-1,3:end-1) + J(:,:,1:end-1);
        end

        if M(i) == 4
            if i ~= 1
                J = Rand >= R;
                R = R + s*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,2:end-1,1:end-2)));
                J = J & Rand < R;
            else
                R = s*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,2:end-1,1:end-2)));
                J = Rand < R;
            end
            Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - J;
            Latt(2:end-1,2:end-1,2:end-2) = Latt(2:end-1,2:end-1,2:end-2) + J(:,:,2:end);
        end

        if M(i) == 5
            if i ~= 1
                J = Rand >= R;
                R = R + u*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(3:end,2:end-1,2:end-1)));
                J = J & Rand < R;
            else
                R = u*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(3:end,2:end-1,2:end-1)));
                J = Rand < R;
            end
            Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - J;
            Latt(3:end-1,2:end-1,2:end-1) = Latt(3:end-1,2:end-1,2:end-1) + J(1:end-1,:,:);
        end

        if M(i) == 6
            if i ~= 1
                J = Rand >= R;
                R = R + v*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(1:end-2,2:end-1,2:end-1)));
                J = J & Rand < R;
            else
                R = v*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(1:end-2,2:end-1,2:end-1)));
                J = Rand < R;
            end
            Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - J;
            Latt(2:end-2,2:end-1,2:end-1) = Latt(2:end-2,2:end-1,2:end-1) + J(2:end,:,:);
        end
    end

    % Deposition

    if Deposition
        SubsTimeCount = SubsTimeCount + 1;

        J0 = Rand(end-SubstrtHeight,SubsDim2,SubsDim3) >= R6(end-SubstrtHeight,SubsDim2,SubsDim3) & Latt(end-SubstrtHeight-1,SubsDim2+1,SubsDim3+1) >= drho/2;
        JD = rand(1,length(SubsDim2),length(SubsDim3)) >= p+q+r+s+u & J0;

        Latt(end-SubstrtHeight-1,SubsDim2+1,SubsDim3+1) = Latt(end-SubstrtHeight-1,SubsDim2+1,SubsDim3+1) - JD*drho;

        Substrate = Substrate + JD*drho;
        
    end

    % Simulaton

    if Simulation && floor(FlowTimeCount/Speed) == FlowTimeCount/Speed

%       TopViewFlowH.ColorData = rot90(squeeze(sum(Latt(2:end-1,2:end-1,2:end-1),1)),1);
%         figure(1)        
%         spy(rot90(squeeze(sum(Latt(2:end-1,2:end-1,2:end-1),1)),1));
%         SidViewFlowH.ColorData = sum(Latt(2:end-1,2:end-1,2:end-1),3);
%         SubstrateH.ColorData = rot90(squeeze(Substrate));
% %       TopViewFlowH.Title = ["Top View","(Time "+ num2str(FlowTimeCount)+")"];
%         SidViewFlowH.Title = ["Side View","(Time "+ num2str(FlowTimeCount)+")"];
%         SubstrateH.Title = ["Deposition","(Time " + num2str(SubsTimeCount)+")"];
        [x,y,z]=ind2sub(size(Latt(2:end-1,2:end-1,2:end-1)),find(Latt(2:end-1,2:end-1,2:end-1)));
        figure(1)        
        set(PlotLatt,'XData',x,'YData',y,'ZData',z)
        drawnow;
    end

    if floor((t/Time)*100) == (t/Time)*100
        clc;
        fprintf('Simulating... (%d%%)\n', (t/Time)*100);
    end

end
toc

% pbaspect
