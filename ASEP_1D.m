%% Simulation Settings

LattS = 100;
[Alfa , Beta] = deal(0.4,0.2);
XBias = 1;     % for TASEP assign 1
Time = 1000;

%% Initiating variables

Randomnes = 1; % Degree of randomnes in bulk
N = 50000;      % numer of times same simulation is done (for ensembled sum)
[p, q] = deal( 1/2+(XBias/2) , 1/2-(XBias/2) );
ExactLattice = zeros(1,LattS);

%% Simulation for Exact Solution
clc
tic
%%
for i = 1:N

    Latt = false(1,LattS);
    % Current [ J(1,:) => towards right, J(2,:) => towars left ] 
    J = false(2,LattS);

    M = randperm(3);

    for t = 1:Time

        for j = 1:3

            %Exclusion
            if M(j) == 1
                % to remove randomnes in remove line 1 2 3
                for k = 1:Randomnes+1                                   % 1
                    R = rand(1,LattS);
                    R(Latt == false) = nan;
                    J(1,1:end-1) = R(1:end-1)<p & ~Latt(2:end);
                    J(2,2:end) = R(2:end)>=p & R(2:end)<p+q & ~Latt(1:end-1);
                    J(:,randperm(LattS,LattS/(Randomnes+1))) = false;   % 2
                    Latt(circshift(J(1,:),1,2)) = true;
                    Latt(circshift(J(2,:),-1,2)) = true;
                    Latt(J(1,:)|J(2,:)) = false;
                end                                                     % 3

            end

            % Injection
            if M(j) == 2
                Latt(1) = Latt(1) | rand < Alfa;
            end

            % Ejection
            if M(j) == 3
                Latt(end) = Latt(end) & ~(rand < Beta);
            end
        end
    end

    if floor((i/N)*100) == (i/N)*100
        clc; fprintf('Simulating... (%d%%)\n', (i/N)*100);
    end

    ExactLattice = ExactLattice + Latt;
end
toc

figure(1)
plot(1:length(ExactLattice),ExactLattice/N);
xlim([0.5 LattS]); ylim([0 1]); grid on;
xlabel("Position ($i$)",'Interpreter','latex',"FontSize",16);
ylabel('$\rho_i$','Interpreter','latex',"FontSize",16)
title('Maximum Current Plot','($\alpha = 0.8,\beta = 0.8$)','Interpreter','latex','FontSize',16)


disp("Exact Soluton ploted")