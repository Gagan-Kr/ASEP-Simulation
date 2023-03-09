%% Lattice Config
LattS = [ 1 100 100 ];   %     for 2D simulation [ 1 Column Row ]            %*

Latt = false(LattS);
J = false(cat(2,LattS,6));

FlowTimeCount = 0;

%% Substrate

[SubsDim2, SubsDim3] =  deal( 21:80 , 1:LattS(3)  );                        %*

Substrate = zeros(1,LattS(2),LattS(3));
Substrate(1,LattS(2),LattS(3)) = 0.001;
LatticeSub = ones(1,length(SubsDim2)+2,length(SubsDim3)+2);
SubsTimeCount = 0;

% f2 = figure(2);
% SubstrateH = heatmap(rot90(squeeze(Substrate),1),"Colormap",turbo);
% [XLabels,YLabels] = deal(string(1:LattS(2)),string(1:LattS(3)));
% XLabels(mod(1:LattS(2),5) ~= 0) = " ";
% YLabels(mod(1:LattS(3),5) ~= 0) = " ";
% set(SubstrateH,"XDisplayLabels",XLabels,"YDisplayLabels",YLabels)
% grid off; caxis([1 100]);

SubstrateEns = rot90(squeeze(Substrate));
% SubstrateEns20000 = rot90(squeeze(Substrate));
% SubstrateEns30000 = rot90(squeeze(Substrate));
% SubstrateEns40000 = rot90(squeeze(Substrate));
% SubstrateEns50000 = rot90(squeeze(Substrate));
% SubstrateEns60000 = rot90(squeeze(Substrate));
% SubstrateEns70000 = rot90(squeeze(Substrate));
% SubstrateEns80000 = rot90(squeeze(Substrate));
% SubstrateEns90000 = rot90(squeeze(Substrate));
% SubstrateEns100000 = rot90(squeeze(Substrate));
% SubstrateEns120000 = rot90(squeeze(Substrate));
% SubstrateEns140000 = rot90(squeeze(Substrate));
% SubstrateEns160000 = rot90(squeeze(Substrate));
% SubstrateEns180000 = rot90(squeeze(Substrate));
% SubstrateEns200000 = rot90(squeeze(Substrate));
% 
% Substrate10000 = rot90(squeeze(Substrate));
% Substrate20000 = rot90(squeeze(Substrate));
% Substrate30000 = rot90(squeeze(Substrate));
% Substrate40000 = rot90(squeeze(Substrate));
% Substrate50000 = rot90(squeeze(Substrate));
% Substrate60000 = rot90(squeeze(Substrate));
% Substrate70000 = rot90(squeeze(Substrate));
% Substrate80000 = rot90(squeeze(Substrate));
% Substrate90000 = rot90(squeeze(Substrate));
% Substrate100000 = rot90(squeeze(Substrate));
% Substrate120000 = rot90(squeeze(Substrate));
% Substrate140000 = rot90(squeeze(Substrate));
% Substrate160000 = rot90(squeeze(Substrate));
% Substrate180000 = rot90(squeeze(Substrate));
% Substrate200000 = rot90(squeeze(Substrate));

% figure(5)
% SubstrateHPlot = plot(1:LattS(2),zeros(1,LattS(2)));
% grid on;
% 
% figure(6)
% SubstrateLPlot = plot(1:LattS(3),zeros(1,LattS(3)));
% grid on;

%% Injection & Ejection window

[IL3, IL1, EL3, EL1] = deal( 46:55 , 1 , [] , [] );                         %*
[IR3, IR1, ER3, ER1] = deal( [] , [] , 1:LattS(3) , 1 );                    %*

%% Simukation Settngs

[ Alfa , Beta ] = deal( 1 , 1 );                                            %*

[ XBias , YBias]  = deal( 0.4 , 0 );                                        %*

Deposition = true;                                                          %*
MaxGrowth = 100; MaxDepRate = 0.03;                                         %*

DiffusionRateSubstrt = 0.005;

Simulation = false; Speed = 1000;                                            %*
Time = 100000;                                                               %*

[p, q] = deal( 1/4+(XBias/4) , 1/4-(XBias/4) );
[r, s] = deal( 1/4+(YBias/4) , 1/4-(YBias/4) );

[ps,qs] = deal(p*DiffusionRateSubstrt,q*DiffusionRateSubstrt);
[rs,ss] = deal(r*DiffusionRateSubstrt,s*DiffusionRateSubstrt);

% Run Simulation
TotRuns = 10;
% figure(1);

tic
for N = 1:TotRuns

    Latt = false(LattS);
    J = false(cat(2,LattS,6));
    FlowTimeCount = 0;

    Substrate = zeros(1,LattS(2),LattS(3));
    Substrate(1,LattS(2),LattS(3)) = 0.001;
    SubsTimeCount = 0;

    for t = 1:Time
        FlowTimeCount = FlowTimeCount + 1;

        % Injection & Ejection

%        Latt(IR1,end,IR3) = Latt(IR1,end,IR3) | rand(length(IR1),1,length(IR3)) < Alfa;
        Latt(IL1,1,IL3) = Latt(IL1,1,IL3) | rand(length(IL1),1,length(IL3)) < Alfa;

        Latt(ER1,end,ER3) = Latt(ER1,end,ER3) & rand(length(ER1),1,length(ER3)) > Beta;
%        Latt(EL1,1,EL3) = Latt(EL1,1,EL3) & rand(length(EL1),1,length(EL3)) > Beta;

        % Deposition

        if Deposition
            SubsTimeCount = SubsTimeCount + 1;

            % Local Deposition rate
            DepRate = MaxDepRate*((MaxGrowth-Substrate(1,SubsDim2,SubsDim3))/MaxGrowth);
            RandDep = rand(1,length(SubsDim2), length(SubsDim3)) < DepRate;
            JD = Latt(1,SubsDim2, SubsDim3) & RandDep;
            Substrate(1,SubsDim2,SubsDim3) = Substrate(1,SubsDim2,SubsDim3) + JD;
            Latt(1,SubsDim2, SubsDim3) = Latt(1,SubsDim2, SubsDim3) & ~ RandDep;

            % Diffusion in Substrate
            LatticeSub(1,2:end-1,2:end-1) = Substrate(1,SubsDim2,SubsDim3)/MaxGrowth;
            R1 = (ps*LatticeSub(1,2:end-1,2:end-1)).*(1-LatticeSub(1,3:end,2:end-1));
            R2 = R1 + ((qs*LatticeSub(1,2:end-1,2:end-1)).*(1-LatticeSub(1,1:end-2,2:end-1)));
            R3 = R2 + ((rs*LatticeSub(1,2:end-1,2:end-1)).*(1-LatticeSub(1,2:end-1,1:end-2)));
            R4 = R3 + ((ss*LatticeSub(1,2:end-1,2:end-1)).*(1-LatticeSub(1,2:end-1,3:end)));

            Rand = rand(1,length(SubsDim2),length(SubsDim3));

            JR = Rand < R1;
            JL = Rand >= R1 & Rand < R2;
            JU = Rand >= R2 & Rand < R3;
            JD = Rand >= R3 & Rand < R4;

            Substrate(1,SubsDim2,SubsDim3) = Substrate(1,SubsDim2,SubsDim3) - (JR+JL+JU+JD);
            Substrate(1,SubsDim2(2):SubsDim2(end),SubsDim3) = Substrate(1,SubsDim2(2):SubsDim2(end),SubsDim3) + JR(1,1:end-1,:);
            Substrate(1,SubsDim2(1):SubsDim2(end-1),SubsDim3) = Substrate(1,SubsDim2(1):SubsDim2(end-1),SubsDim3) + JL(1,2:end,:);
            Substrate(1,SubsDim2,SubsDim3(1):SubsDim3(end-1)) = Substrate(1,SubsDim2,SubsDim3(1):SubsDim3(end-1)) + JU(1,:,2:end);
            Substrate(1,SubsDim2,SubsDim3(2):SubsDim3(end)) = Substrate(1,SubsDim2,SubsDim3(2):SubsDim3(end)) + JD(1,:,1:end-1);

            % Injection from substrate
            Substrate(1,SubsDim2(end),SubsDim3) = Substrate(1,SubsDim2(end),SubsDim3) - (rand(size(Substrate(1,SubsDim2(end),SubsDim3)))<Beta.*(Substrate(1,SubsDim2(end),SubsDim3)/MaxGrowth));

        end

% %         Asymetric Exclusion Process
% % 
% %         R = rand(LattS);
% %         R(Latt == false) = nan;
% %         M = randperm(4);
% % 
% %         for i = 1:4
% % 
% %             if M(i) == 1   % Right hop
% %                 J(:,1:end-1,:,1) = R(:,1:end-1,:)<p & ~Latt(:,2:end,:);
% %                 Latt(circshift(J(:,:,:,1),1,2)) = true;
% %             end
% % 
% %             if M(i) == 2   % Left hop
% %                 J(:,2:end,:,2) = R(:,2:end,:)>=p & R(:,2:end,:)<p+q & ~Latt(:,1:end-1,:);
% %                 Latt(circshift(J(:,:,:,2),-1,2)) = true;
% %             end
% % 
% %             if M(i) == 3   % Up hop
% %                 J(:,:,1:end-1,3) = (R(:,:,1:end-1)>=p+q & ~Latt(:,:,2:end) & R(:,:,1:end-1)<p+q+r);
% %                 Latt(circshift(J(:,:,:,3),1,3)) = true;
% %             end
% % 
% %             if M(i) == 4   % Down hop
% %                 J(:,:,2:end,4) = R(:,:,2:end)>=p+q+r & ~Latt(:,:,1:end-1) & R(:,:,2:end)<p+q+r+s;
% %                 Latt(circshift(J(:,:,:,4),-1,3)) = true;
% %             end
% %         end
% % 
% %         Latt(J(:,:,:,1)|J(:,:,:,2)|J(:,:,:,3)|J(:,:,:,4)) = false;

        % Optimisation Test

        % Time evolution issue, but much faster (5 times) than previous algorithm

        R = rand;

        if R < p                            % Move right
            J(:,1:end-1,:,1) = Latt(:,1:end-1,:) & ~Latt(:,2:end,:);
            Latt(circshift(J(:,:,:,1),1,2)) = true;
            Latt(J(:,:,:,1)) = false;
        elseif R > p && R < p+q             % Move down
            J(:,2:end,:,2) = Latt(:,2:end,:) & ~Latt(:,1:end-1,:);
            Latt(circshift(J(:,:,:,2),-1,2)) = true;
            Latt(J(:,:,:,2)) = false;
            
        elseif R > p+q && R < p+q+r         % Move up
            J(:,:,1:end-1,3) = Latt(:,:,1:end-1) & ~Latt(:,:,2:end);
            Latt(circshift(J(:,:,:,3),1,3)) = true;
            Latt(J(:,:,:,3)) = false;
        elseif R > p+q+r                    % Move down
            J(:,:,2:end,4) = Latt(:,:,2:end) & ~Latt(:,:,1:end-1);
            Latt(circshift(J(:,:,:,4),-1,3)) = true;
            Latt(J(:,:,:,4)) = false;
        end


        % Simulaton

        if Simulation && floor(FlowTimeCount/Speed) == FlowTimeCount/Speed

            spy(rot90(squeeze(Latt),1)); grid on;
            title("Particle flow","(Time "+num2str(FlowTimeCount)+")");
            SubstrateData = rot90(squeeze(Substrate));
            set(SubstrateH,'ColorData',SubstrateData);
            SubstrateH.Title = ["Deposition","(Time "+num2str(SubsTimeCount)+")"];

            Substrate2dH = sum(SubstrateData((LattS(3)/2)-1:LattS(3)/2+2,:),1)/4;
            [ ~ , Ind] = max(Substrate2dH);
            if Ind < 2
                Ind = 3;
                Substrate2dV = sum(SubstrateData(:,Ind-1:Ind+1),2)/3;
            else
                Substrate2dV = sum(SubstrateData(:,Ind-1:Ind+1),2)/3;
            end

            set(SubstrateHPlot,'YData',Substrate2dH);
            set(SubstrateLPlot,'YData',Substrate2dV);

            drawnow;
        end

%         if t == 10000
%             SubstrateEns10000 = SubstrateEns10000 + rot90(squeeze(Substrate));
%             Substrate10000 = rot90(squeeze(Substrate));
%         elseif t == 20000
%             SubstrateEns20000 = SubstrateEns20000 + rot90(squeeze(Substrate));
%             Substrate20000 = rot90(squeeze(Substrate));
%         elseif t == 30000
%             SubstrateEns30000 = SubstrateEns30000 + rot90(squeeze(Substrate));
%             Substrate30000 = rot90(squeeze(Substrate));
%         elseif t == 40000
%             SubstrateEns40000 = SubstrateEns40000 + rot90(squeeze(Substrate));
%             Substrate40000 = rot90(squeeze(Substrate));
%         elseif t == 50000
%             SubstrateEns50000 = SubstrateEns50000 + rot90(squeeze(Substrate));
%             Substrate50000 = rot90(squeeze(Substrate));
%         elseif t == 60000
%             SubstrateEns60000 = SubstrateEns60000 + rot90(squeeze(Substrate));
%             Substrate60000 = rot90(squeeze(Substrate));
%         elseif t == 70000
%             SubstrateEns70000 = SubstrateEns70000 + rot90(squeeze(Substrate));
%             Substrate70000 = rot90(squeeze(Substrate));
%         elseif t == 80000
%             SubstrateEns80000 = SubstrateEns80000 + rot90(squeeze(Substrate));
%             Substrate80000 = rot90(squeeze(Substrate));
%         elseif t == 90000
%             SubstrateEns90000 = SubstrateEns90000 + rot90(squeeze(Substrate));
%             Substrate90000 = rot90(squeeze(Substrate));
%         elseif t == 100000
%             SubstrateEns100000 = SubstrateEns100000 + rot90(squeeze(Substrate));
%             Substrate100000 = rot90(squeeze(Substrate));
%         elseif t == 120000
%             SubstrateEns120000 = SubstrateEns120000 + rot90(squeeze(Substrate));
%             Substrate120000 = rot90(squeeze(Substrate));
%         elseif t == 140000
%             SubstrateEns140000 = SubstrateEns140000 + rot90(squeeze(Substrate));
%             Substrate140000 = rot90(squeeze(Substrate));
%         elseif t == 160000
%             SubstrateEns160000 = SubstrateEns160000 + rot90(squeeze(Substrate));
%             Substrate160000 = rot90(squeeze(Substrate));
%         elseif t == 180000
%             SubstrateEns180000 = SubstrateEns180000 + rot90(squeeze(Substrate));
%             Substrate180000 = rot90(squeeze(Substrate));
%         elseif t == 200000
%             SubstrateEns200000 = SubstrateEns200000 + rot90(squeeze(Substrate));
%             Substrate200000 = rot90(squeeze(Substrate));
%         end

        if floor((t/Time)*100) == (t/Time)*100
            clc; fprintf('Simulating... (%d/%d) (%d%%)\n', N, TotRuns, (t/Time)*100);
        end
    end

    SubstrateEns = SubstrateEns + rot90(squeeze(Substrate));
end
toc

SubstrateEnsRB4 = SubstrateEns;
% SubstrateEns10000 = SubstrateEns10000/TotRuns;
% SubstrateEns20000 = SubstrateEns20000/TotRuns;
% SubstrateEns30000 = SubstrateEns30000/TotRuns;
% SubstrateEns40000 = SubstrateEns40000/TotRuns;
% SubstrateEns50000 = SubstrateEns50000/TotRuns;
% SubstrateEns60000 = SubstrateEns60000/TotRuns;
% SubstrateEns70000 = SubstrateEns70000/TotRuns;
% SubstrateEns80000 = SubstrateEns80000/TotRuns;
% SubstrateEns90000 = SubstrateEns90000/TotRuns;
% SubstrateEns100000 = SubstrateEns100000/TotRuns;
% SubstrateEns120000 = SubstrateEns120000/TotRuns;
% SubstrateEns140000 = SubstrateEns140000/TotRuns;
% SubstrateEns160000 = SubstrateEns160000/TotRuns;
% SubstrateEns180000 = SubstrateEns180000/TotRuns;
% SubstrateEns200000 = SubstrateEns200000/TotRuns;







% %% Lattice Config
% 
% LattSize = [ 30 50 30 ];    %  [ Height Length Width ]
% SubLattSize = 125;
% 
% FlowTimeCount = 0;
% drho = 1/SubLattSize;
% Latt = zeros(LattSize+2);
% Latt(1,:,:) = 1; Latt(end,:,:) = 1;
% Latt(:,1,:) = 1; Latt(:,end,:) = 1;
% Latt(:,:,1) = 1; Latt(:,:,end) = 1;
% 
% %% Injection & Ejection window
% 
% [InjL3, InjL1, EjeL3, EjeL1] = deal( 11:20 , 13:18 , [] , [] );
% [InjR3, InjR1, EjeR3, EjeR1] = deal( [] , [] , 1:30 , 1:30 );
% [InjU3, InjU2, EjeU3, EjeU2] = deal( [] , [] , [] , [] );
% [InjD3, InjD2, EjeD3, EjeD2] = deal( [] , [] , [] , [] );
% 
% figure(1)
% TopViewFlowH = heatmap(rot90(squeeze(sum(Latt(2:end-1,2:end-1,2:end-1),1)),1),"Colormap",jet);
% TopViewFlowH.Title = ["Top View","(Time "+ num2str(FlowTimeCount)+")"];
% grid off;
% 
% figure(2)
% SidViewFlowH = heatmap(sum(Latt(2:end-1,2:end-1,2:end-1),3),"Colormap",jet);
% SidViewFlowH.Title = ["Side View","(Time "+ num2str(FlowTimeCount)+")"];
% grid off;
% 
% %% Substrate
% 
% [SubsDim2, SubsDim3] =  deal( 1:LattSize(2) , 1:LattSize(3) );
% 
% Substrate = zeros(1,LattSize(2),LattSize(3));
% Substrate(1,LattSize(2),LattSize(3)) = 0.001;
% SubsTimeCount = 0;
% 
% figure(3)
% SubstrateH = heatmap(rot90(squeeze(Substrate),1),"Colormap",turbo);
% grid off;
% 
% %% Simukation Settngs
% 
% Alfa =  1 ;
% Beta = 0.7;
% 
% XBias = 0.2;
% YBias =  0 ;
% ZBias =  0 ;
% 
% [p, q] = deal( 1/6+(XBias/6) , 1/6-(XBias/6) );
% [r, s] = deal( 1/6+(YBias/6) , 1/6-(YBias/6) );
% [u, v] = deal( 1/6+(ZBias/6) , 1/6-(ZBias/6) );
% 
% Deposition = true; RxnRate = 1;
% 
% Simulation = true; Speed = 100;
% 
% Time = 10000;
% 
% %% Run Simulation
% tic 
% for t = 1:Time
% 
%     FlowTimeCount = FlowTimeCount + 1;
% 
%     % Injection & Ejection
% 
%     AlfaRateR = Alfa*(1-Latt(InjR1+1,end-1,InjR3+1));
%     RandInjR = rand(size(AlfaRateR));
%     Latt(InjR1+1,end-1,InjR3+1) = Latt(InjR1+1,end-1,InjR3+1) + drho*(RandInjR<AlfaRateR);
% 
%     AlfaRateD = Alfa*(1-Latt(end-1,InjD2+1,InjD3+1));
%     RandInjD = rand(size(AlfaRateD));
%     Latt(end-1,InjD2+1,InjD3+1) = Latt(end-1,InjD2+1,InjD3+1) + drho*(RandInjD<AlfaRateD);
% 
%     AlfaRateL = Alfa*(1-Latt(InjL1+1,2,InjL3+1));
%     RandInjL = rand(size(AlfaRateL));
%     Latt(InjL1+1,2,InjL3+1) = Latt(InjL1+1,2,InjL3+1) + drho*(RandInjL<AlfaRateL);
% 
%     AlfaRateU = Alfa*(1-Latt(2,InjU2+1,InjU3+1));
%     RandInjU = rand(size(AlfaRateU));
%     Latt(2,InjU2+1,InjU3+1) = Latt(2,InjU2+1,InjU3+1) + drho*(RandInjU<AlfaRateU);
% 
%     BetaRateR = Beta*Latt(EjeR1+1,end-1,EjeR3+1);
%     RandEjeR = rand(size(BetaRateR));
%     Latt(EjeR1+1,end-1,EjeR3+1) = Latt(EjeR1+1,end-1,EjeR3+1) - drho*(RandEjeR<BetaRateR);
% 
%     BetaRateD = Beta*Latt(end-1,EjeD2+1,EjeD3+1);
%     RandEjeD = rand(size(BetaRateD));
%     Latt(end-1,EjeD2+1,EjeD3+1) = Latt(end-1,EjeD2+1,EjeD3+1) - drho*(RandEjeD<BetaRateD);
% 
%     BetaRateL = Beta*Latt(EjeL1+1,2,EjeL3+1);
%     RandEjeL = rand(size(BetaRateL));
%     Latt(EjeL1+1,2,EjeL3+1) = Latt(EjeL1+1,2,EjeL3+1) - drho*(RandEjeL<BetaRateL);
% 
%     BetaRateU = Beta*Latt(2,EjeU2+1,EjeL3+1);
%     RandEjeU = rand(size(BetaRateU));
%     Latt(2,EjeU2+1,EjeU3+1) = Latt(2,EjeU2+1,EjeU3+1) - drho*(RandEjeU<BetaRateU);
% 
%     % Asymetric Exclusion Process
% 
%     Rand = rand(LattSize);
%     M = randperm(6);
% 
%     for i = 1:6
% 
%         if M(i) == 1
%             if i ~= 1
%                 J = Rand >= R;
%                 R = R + p*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,3:end,2:end-1)));
%                 J = J & Rand < R;
%             else
%                 R = p*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,3:end,2:end-1)));
%                 J = Rand < R;
%             end
%             Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - J*drho;
%             Latt(2:end-1,3:end-1,2:end-1) = Latt(2:end-1,3:end-1,2:end-1) + J(:,1:end-1,:)*drho;
%         end
% 
%         if M(i) == 2
%             if i ~= 1
%                 J = Rand >= R;
%                 R = R + q*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,1:end-2,2:end-1)));
%                 J = J & Rand < R;
%             else
%                 R = q*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,1:end-2,2:end-1)));
%                 J = Rand < R;
%             end
%             Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - J*drho;
%             Latt(2:end-1,2:end-2,2:end-1) = Latt(2:end-1,2:end-2,2:end-1) + J(:,2:end,:)*drho;
%         end
%  
%         if M(i) == 3
%             if i ~= 1
%                 J = Rand >= R;
%                 R = R + r*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,2:end-1,3:end)));
%                 J = J & Rand < R;
%             else
%                 R = r*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,2:end-1,3:end)));
%                 J = Rand < R;
%             end
%             Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - drho*J;
%             Latt(2:end-1,2:end-1,3:end-1) = Latt(2:end-1,2:end-1,3:end-1) + drho*J(:,:,1:end-1);
%         end
% 
%         if M(i) == 4
%             if i ~= 1
%                 J = Rand >= R;
%                 R = R + s*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,2:end-1,1:end-2)));
%                 J = J & Rand < R;
%             else
%                 R = s*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(2:end-1,2:end-1,1:end-2)));
%                 J = Rand < R;
%             end
%             Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - drho*J;
%             Latt(2:end-1,2:end-1,2:end-2) = Latt(2:end-1,2:end-1,2:end-2) + drho*J(:,:,2:end);
%         end
% 
%         if M(i) == 5
%             if i ~= 1
%                 J = Rand >= R;
%                 R = R + u*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(3:end,2:end-1,2:end-1)));
%                 J = J & Rand < R;
%             else
%                 R = u*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(3:end,2:end-1,2:end-1)));
%                 J = Rand < R;
%             end
%             Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - drho*J;
%             Latt(3:end-1,2:end-1,2:end-1) = Latt(3:end-1,2:end-1,2:end-1) + drho*J(1:end-1,:,:);
%         end
% 
%         if M(i) == 6
%             if i ~= 1
%                 J = Rand >= R;
%                 R = R + v*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(1:end-2,2:end-1,2:end-1)));
%                 J = J & Rand < R;
%             else
%                 R = v*(Latt(2:end-1,2:end-1,2:end-1).*(1-Latt(1:end-2,2:end-1,2:end-1)));
%                 J = Rand < R;
%             end
%             Latt(2:end-1,2:end-1,2:end-1) = Latt(2:end-1,2:end-1,2:end-1) - drho*J;
%             Latt(2:end-2,2:end-1,2:end-1) = Latt(2:end-2,2:end-1,2:end-1) + drho*J(2:end,:,:);
%         end
%     end
% 
%     % Deposition
% 
%     if Deposition && rand < RxnRate
%         SubsTimeCount = SubsTimeCount + 1;
%         J0 = Rand(end,SubsDim2,SubsDim3) >= R(end,SubsDim2,SubsDim3) & Latt(end-1,SubsDim2+1,SubsDim3+1) >= drho/2;
%         JD = rand(1,length(SubsDim2),length(SubsDim3)) >= p+q+r+s+u & J0;
% 
%         Latt(end-1,SubsDim2+1,SubsDim3+1) = Latt(end-1,SubsDim2+1,SubsDim3+1) - JD*drho;
%         Substrate = Substrate + JD*drho;
%     end
% 
%     % Simulaton
% 
%     if Simulation && floor(FlowTimeCount/Speed) == FlowTimeCount/Speed
% 
%         TopViewFlowH.ColorData = rot90(squeeze(sum(Latt(2:end-1,2:end-1,2:end-1),1)),1);
%         SidViewFlowH.ColorData = sum(Latt(2:end-1,2:end-1,2:end-1),3);
%         SubstrateH.ColorData = rot90(squeeze(Substrate));
%         TopViewFlowH.Title = ["Top View","(Time "+ num2str(FlowTimeCount)+")"];
%         SidViewFlowH.Title = ["Side View","(Time "+ num2str(FlowTimeCount)+")"];
%         SubstrateH.Title = ["Deposition","(Time " + num2str(SubsTimeCount)+")"];
%         drawnow;
%     end
% 
%     if floor((t/Time)*100) == (t/Time)*100
%         clc;
%         fprintf('Simulating... (%d%%)\n', (t/Time)*100);
%     end
% end
% toc
% 


%%%% GIF Producer

% Gif code is temporarily deleted
