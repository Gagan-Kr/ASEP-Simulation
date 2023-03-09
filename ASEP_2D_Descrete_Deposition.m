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

f2 = figure(2);
SubstrateH = heatmap(rot90(squeeze(Substrate),1),"Colormap",turbo);
[XLabels,YLabels] = deal(string(1:LattS(2)),string(1:LattS(3)));
XLabels(mod(1:LattS(2),5) ~= 0) = " ";
YLabels(mod(1:LattS(3),5) ~= 0) = " ";
set(SubstrateH,"XDisplayLabels",XLabels,"YDisplayLabels",YLabels)
grid off; caxis([1 100]);

SubstrateEns = rot90(squeeze(Substrate));


figure(3)
SubstrateHPlot = plot(1:LattS(2),zeros(1,LattS(2)));
grid on;

figure(4)
SubstrateLPlot = plot(1:LattS(3),zeros(1,LattS(3)));
grid on;

%% Injection & Ejection window

[IL3, IL1, EL3, EL1] = deal( 46:55 , 1 , [] , [] );                         %*
[IR3, IR1, ER3, ER1] = deal( [] , [] , 1:LattS(3) , 1 );                    %*

%% Simukation Settngs

[ Alfa , Beta ] = deal( 1 , 1 );                                            %*

[ XBias , YBias]  = deal( 0.4 , 0 );                                        %*

Deposition = true;                                                          %*
MaxGrowth = 100; MaxDepRate = 0.03;                                         %*

DiffusionRateSubstrt = 0.05;

Simulation = true; Speed = 100;                                            %*
Time = 10000;                                                               %*

[p, q] = deal( 1/4+(XBias/4) , 1/4-(XBias/4) );
[r, s] = deal( 1/4+(YBias/4) , 1/4-(YBias/4) );

[ps,qs] = deal(p*DiffusionRateSubstrt,q*DiffusionRateSubstrt);
[rs,ss] = deal(r*DiffusionRateSubstrt,s*DiffusionRateSubstrt);

% Run Simulation
TotRuns = 10;
figure(1);

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
        R = rand(LattS);
        R(Latt == false) = nan;
        M = randperm(4);

        for i = 1:4

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
        end

        Latt(J(:,:,:,1)|J(:,:,:,2)|J(:,:,:,3)|J(:,:,:,4)) = false;

        % Optimisation Test

%         R = rand;
% 
%         if R < p                            % Move right
%             J(:,1:end-1,:,1) = Latt(:,1:end-1,:) & ~Latt(:,2:end,:);
%             Latt(circshift(J(:,:,:,1),1,2)) = true;
%             Latt(J(:,:,:,1)) = false;
%         elseif R > p && R < p+q                % Move down
%             J(:,2:end,:,2) = Latt(:,2:end,:) & ~Latt(:,1:end-1,:);
%             Latt(circshift(J(:,:,:,2),-1,2)) = true;
%             Latt(J(:,:,:,2)) = false;
%             
%         elseif R > p+q && R < p+q+r         % Move up
%             J(:,:,1:end-1,3) = Latt(:,:,1:end-1) & ~Latt(:,:,2:end);
%             Latt(circshift(J(:,:,:,3),1,3)) = true;
%             Latt(J(:,:,:,3)) = false;
%         elseif R > p+q+r                    % Move down
%             J(:,:,2:end,4) = Latt(:,:,2:end) & ~Latt(:,:,1:end-1);
%             Latt(circshift(J(:,:,:,4),-1,3)) = true;
%             Latt(J(:,:,:,4)) = false;
%         end

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

        if floor((t/Time)*100) == (t/Time)*100
            clc; fprintf('Simulating... (%d/%d) (%d%%)\n', N, TotRuns, (t/Time)*100);
        end
    end

    SubstrateEns = SubstrateEns + rot90(squeeze(Substrate));
end
toc

SubstrateEnsD05 = SubstrateEns/TotRuns;

clf; heatmap(SubstrateEns(:,21:80)/TotRuns,"Colormap",turbo);
grid off; caxis([1 100]);
