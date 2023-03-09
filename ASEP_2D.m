%% Lattice Config
LatticeSize =     [70 , 70];
SubLatticeSize = 200;

Lattice = zeros(LatticeSize);

LatticeB = zeros(LatticeSize+2);
LatticeB(:,1) = 1; LatticeB(:,end) = 1;
LatticeB(1,:) = 1; LatticeB(end,:) = 1;

Count = 0;
FirstRun = true;
FlowTimeCounter = 0;
drho = 1/SubLatticeSize;


%% Injection & Ejection

[InjR, InjL, InjU, InjD] = deal([],21:50,[],[]);
[EjeR, EjeL, EjeU, EjeD] = deal(1:70,[],[],[]);

%% Substrate

FromRowToRow = 31:40;
FromColToCol = 21:50;

[SR,SC] = deal( FromRowToRow , FromColToCol );

Substrate = zeros(LatticeSize);
Substrate(LatticeSize(1),LatticeSize(2)) = 0.0001;
SubstrateTimeCounter = 0;

%% Simulation Settings

Alfa = 1;
Beta = 1;

RightBias = 0.4;
[p, q, r, s] = deal( (0.5+(RightBias)/2)*(1/2) , (0.5-(RightBias)/2)*(1/2)  , 1/4 , 1/4);

Deposition = true; ReactionRate = 0.7; ThicSatu = 10;
Simulation = true; Speed = 10;

Time = 1000;

warning('off','all')
exist FirstRun;

if ans == 0
    error(['Initialize system configuration settings' ...
        ' by pressing (Initialize/Reset Flow System) button :)'])
end

exist Substrate;

if ans == 0
    error(['Initialize substrate configuration settings by ' ...
        'pressing (Initialize/Reset Substrate) button :)'])
end

exist Alfa;

if ans == 0
    error(['Initialize Simulation Settins by cheking and then ' ...
        'uncheking the Simulation checkbox or vise versa :)'])
end

if FirstRun
    figure(1)
    
    subplot(2,2,1)
    subplot('Position',[0.01,0.34,0.37,0.46]);
    hl = heatmap(Lattice,"Colormap",jet);
    hl.FontColor = 'white';
    caxis([0, 1]); grid off

    subplot(2,2,2)
    subplot('Position',[0.50,0.34,0.38,0.46]);
    hs = heatmap(Substrate,'Colormap',turbo);
    hs.FontColor = 'white';
    grid off

    f = gcf;
    f.Color = [0 0.03 0.1];

    SavedLattice = zeros(LatticeSize(1),LatticeSize(2),3);
    RecordFlowTimeCounter = zeros(1,3);
    Observer = zeros(length(SR),length(SC));

    FirstRun = false;
end
 
ReactionTimeStep = SubLatticeSize - (ReactionRate*SubLatticeSize);

tic

for t = 1:Time*SubLatticeSize

    FlowTimeCounter = FlowTimeCounter + 1;
    
    LatticeB(2:end-1,2:end-1) = Lattice;

    R1 = (p*LatticeB(2:end-1,2:end-1)).*(1-LatticeB(2:end-1,3:end));
    R2 = R1 + ((q*LatticeB(2:end-1,2:end-1)).*(1-LatticeB(2:end-1,1:end-2)));
    R3 = R2 + ((r*LatticeB(2:end-1,2:end-1)).*(1-LatticeB(1:end-2,2:end-1)));
    R4 = R3 + ((s*LatticeB(2:end-1,2:end-1)).*(1-LatticeB(3:end,2:end-1)));

    Rand = rand(LatticeSize);

    JR = Rand < R1;
    JL = Rand >= R1 & Rand < R2;
    JU = Rand >= R2 & Rand < R3;
    JD = Rand >= R3 & Rand < R4;
    JO = Rand >= R4 & Lattice > 0;

    if Deposition

        StayProb = max(R4,[],"all") - R4;

        SubstrateTimeCounter = SubstrateTimeCounter + 1;
%        Observer = JO(SR,SC).*(Observer + 1);

        RandDep = rand(length(SR),length(SC));
%        S = Observer >= ReactionTimeStep & RandDep < (ThicSatu*exp(-Substrate(SR,SC)));
        S = RandDep < Lattice(SR,SC);
        Substrate(SR,SC) = Substrate(SR,SC) + (drho*(S));
        Lattice(SR,SC) = Lattice(SR,SC) - (drho*(S));

%        Observer(S) = 0;

    end    

    Lattice = Injection(Lattice,drho,Alfa,InjR,InjL,InjU,InjD);
    Lattice = Ejection( Lattice,drho,Beta,EjeR,EjeL,EjeU,EjeD);

    Lattice = Lattice - (JR+JL+JU+JD)*drho;
    Lattice(:,2:end) = Lattice(:,2:end) + (drho*JR(:,1:end-1));
    Lattice(:,1:end-1) = Lattice(:,1:end-1) + drho*JL(:,2:end);
    Lattice(1:end-1,:) = Lattice(1:end-1,:) + drho*JU(2:end,:);
    Lattice(2:end,:)   = Lattice(2:end,:) + drho*JD(1:end-1,:);

    SpeedVariable = (FlowTimeCounter/SubLatticeSize)/Speed;

    if t == 1
            fprintf(['Simulating...                               ' ...
                '                                              /- 100 %%'])
    end

    if Simulation && floor(SpeedVariable) == SpeedVariable
        hl.ColorData = Lattice;
        hl.Title = ['Time ',num2str((FlowTimeCounter/SubLatticeSize))];
        hs.Title = ['Deposition Time ',num2str((SubstrateTimeCounter/SubLatticeSize))];
        hs.ColorData = Substrate;
        drawnow;

    end

    CompletionVariable = (t/(Time*SubLatticeSize))*100;
     
    if floor(CompletionVariable) == CompletionVariable
            fprintf('=')
    end
end


function Lattice = Injection(Lattice,drho,Alfa,InjRight,InjLeft,InjUp,InjDown)

    if ~isempty(InjRight)
        AlfaRateMatRight = Alfa*(1-Lattice(InjRight,end));
        RandInjR = rand(size(AlfaRateMatRight));
        Lattice(InjRight,end) = Lattice(InjRight,end) + drho*(RandInjR<AlfaRateMatRight);
    end

    if ~isempty(InjLeft)
        AlfaRateMatLeft = Alfa*(1-Lattice(InjLeft,1));
        RandInjL = rand(size(AlfaRateMatLeft));
        Lattice(InjLeft,1) = Lattice(InjLeft,1) + drho*(RandInjL<AlfaRateMatLeft);
    end

    if ~isempty(InjUp)
        AlfaRateMatUp = Alfa*(1-Lattice(1,InjUp));
        RandInjU = rand(size(AlfaRateMatUp));
        Lattice(1,InjUp) = Lattice(1,InjUp) + drho*(RandInjU<AlfaRateMatUp);
    end

    if ~isempty(InjDown)
        AlfaRateMatDown = Alfa*(1-Lattice(end,InjDown));
        RandInjD = rand(size(AlfaRateMatDown));
        Lattice(end,InjDown) = Lattice(end,InjDown) + drho*(RandInjD<AlfaRateMatDown);
    end
end

function Lattice = Ejection(Lattice,drho,Beta,EjeRight,EjeLeft,EjeUp,EjeDown)

    if ~isempty(EjeRight)
        BetaRateMatRight = Beta.*Lattice(EjeRight,end);
        RandEjeR = rand(size(BetaRateMatRight));
        Lattice(EjeRight,end) = Lattice(EjeRight,end) - (drho*(RandEjeR<BetaRateMatRight));
    end

    if ~isempty(EjeLeft)
        BetaRateMatLeft = Beta*Lattice(EjeLeft,1);
        RandEjeL = rand(size(BetaRateMatLeft));
        Lattice(EjeLeft,1) = Lattice(EjeLeft,1) - (drho*(RandEjeL<BetaRateMatLeft));
    end

    if ~isempty(EjeUp)
        BetaRateMatUp = Beta*Lattice(1,EjeUp);
        RandEjeU = rand(size(BetaRateMatUp));
        Lattice(1,EjeUp) = Lattice(1,EjeUp) - (drho*(RandEjeU<BetaRateMatUp));
    end

    if ~isempty(EjeDown)
        BetaRateMatDown = Beta.*Lattice(end,EjeDown);
        RandEjeD = rand(size(BetaRateMatDown));
        Lattice(end,EjeDown) = Lattice(end,EjeDown) - (drho*(RandEjeD<BetaRateMatDown));
    end
end
