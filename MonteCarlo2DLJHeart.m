

function [finalU,finalVirial,finalPressure,finalConfiguration,finalDistances,moveCount] = ...
    MonteCarlo2DLJHeart(N,T,rho,~,maxdr,initialConfig,~...
    ,initialDistances,initialU,varargin)

%% Monte-Carlo in NVT ensemble for Lennard-Jonse potantioal in 2D %%

% calculate the energy after 'Nstep' Monte-Carlo runs (Metropolis algorithm)
% for a system of 'N' particles in 2D interacting via Lennard-Jonse pair
% potantial. periodic boundary conditions (PBC) are apllied.

% inputs:
% N - number of particles
% T - reduced temperature
% Nsteps - number of steps
% maxdr - maximum particle displacement
% initialConfig - initial configuration of particles (2 by N matrix)
% initialU - initial energy of the configuration
% rCutoff - the cutoff distance for the energy
% optional: the 'm' in the intetraction (U = 4((1/r)^12 - (1/r)^-m)),
%           default is 6.


% outputs:
% finalU - the energy in each step (1 by Nstep matrix)
% finalConfigurations - the coordinates of all particles in each step (2 by N
%                     by Nstep matrix)
% finalDistances - the pair distances of all particles in each step (N by N
%                by Nstep matrix)
% moveCount - counts accepted moves

% the potantial of Lennard-Jonse in reduced units:
%   U = 4*[(1/r)^12 - (1/r)^6] 
%   (you can change the 6 with the optional input 'm')

% the potantial of Lennard-Jonse in non-reduced units:
%   U = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]

% reduced units:
% T(reduced) = kT/epsilon | r(reduced) = r/sigma | U(reduced) = U/epsilon

% usage examples: 

% create an initial configuration of 100 particles in a 100*100 box:
N = 100; 
L = 100;
initialConfig = L*rand(2,N);
rho = N/L^2;


% choose simulation parameters:
 T = 0.727; Nsteps = 1; maxdr = 1; rCutoff = 2.5;


% calculate the pair distances: 
 initialDistances = sqrt(bsxfun(@(x1,x2) (x1-x2).^2 ,...
       initialConfig(1,:),initialConfig(1,:)')...
       +bsxfun(@(x3,x4) (x4-x3).^2 ,...
       initialConfig(2,:),initialConfig(2,:)'));

   
% calculate the initial energy (up to the cutoff):
 d = initialDistances(and(initialDistances < rCutoff,initialDistances > 0));
 initialU = 4*sum(sum(d.^(-12)-d.^(-6)));


% calculate energy
% [finalU,~,~,finalConfiguration,finalDistances,moveCount] = ...
 %    MonteCarlo2DLJHeart(N,T,rho,Nsteps,maxdr,initialConfig,rCutoff...
  %   ,initialDistances,initialU)

% also calculate the pressure.
% first calculate the initial virial:
%virial = -(2*rho/N)*(sum(sum(6*d^(-6) - 12*d^(-12))));

% now calculate the energy and pressure using monte carlo
 %[finalU,finalVirial,finalPressure,finalConfiguration,finalDistances,moveCount] = ...
  %   MonteCarlo2DLJHeart(N,T,rho,Nsteps,maxdr,initialConfig,rCutoff...
   %  ,initialDistances,initialU,'virial',virial);
 


% parse input parameters
p = inputParser();
addOptional(p, 'virial', []);
addOptional(p, 'm', 6);
parse(p, varargin{:});
Results = p.Results;
virial = Results.virial;
m = Results.m;

% initiate variables
dist = initialDistances ;
particlesPosition = initialConfig;
U = initialU;
if ~isempty(virial)
    V = virial;
end
L = sqrt(N/rho); % board length in reduced units
moveCount = 0;
movedParticle = 0;

for step = 1:Nsteps
    
        % choose particle to move 
        movedParticle = movedParticle + 1;
        if movedParticle == N + 1
            movedParticle = 1;
        end
               
                            
        % choose displacement:
        displacex = maxdr*rand - (maxdr/2);
        displacey = maxdr*rand - (maxdr/2);
        displace = sqrt(displacex^2 + displacey^2);

        % move particle
        newParticlesPosition = movePBC(particlesPosition,movedParticle,...
            displacex,displacey,L);

        % calculate new distances
        newDist = reCalcDist(dist,movedParticle,...
            newParticlesPosition,N,L,[]);
        
        % calculate the change in energy
        dU = Uchange(movedParticle,dist,newDist,N,rCutoff,m);
           
       
        % calculate the change in the virial 
        if ~isempty(virial)
            dV = Vchange(movedParticle,dist,newDist,N,rCutoff,rho,m);
        end
        
        % if dU < 0 eccept move
        %(if (1/T)*dU > 75, we are sure the move
        % will not be excepted, so we don't calculate exp(1/T)*dU to 
        % save calculation time)
            
        if (1/T)*dU < 75
            if dU < 0  
                U = U + dU;
                if ~isempty(virial)
                    V = V + dV;
                end
                dist = newDist;
                particlesPosition = newParticlesPosition;
                moveCount = moveCount + 1; 
                        
            else
                %% otherwise,
                % keep the new state with a probability corresponding to the
                % Boltzmann factor. if the new state is rejected, recount the
                % old configuration. 

                if rand < exp(-(1/T)*dU)
                    U = U + dU;
                    
                    
                    if ~isempty(virial)
                        V = V + dV;
                    end
                    dist = newDist;
                    particlesPosition = newParticlesPosition;
                    moveCount = moveCount + 1;
             
                    
                end
            end
        end
end

finalU = U;
    
  

if ~isempty(virial)
    finalVirial = V;
    finalPressure = T*rho + V;
else
    finalVirial = [];
    finalPressure = [];
end
finalConfiguration = particlesPosition;
finalDistances = dist;



%% functions used in main code %%

        function newparticlesPosition = movePBC(particlesPosition,...
                movedParticle,displacex,displacey,L)
            
                % if the particle gets out of the board, apply PBC
                newparticlesPosition = particlesPosition;
                x = particlesPosition(1,movedParticle)+ displacex;
                y = particlesPosition(2,movedParticle)+ displacey;
                
                if x > L/2
                    x = x - L;
                end
                
                if x < (-L/2) 
                    x =  x + L;
                end
                
                if y > L/2 
                    y = y - L;
                end
                
                if y < (-L/2) 
                    y = y + L;
                end
                
                newparticlesPosition(1,movedParticle) = x;
                newparticlesPosition(2,movedParticle) = y;
        end
    
        function newdist = reCalcDist(dist,movedParticles,...
                newParticlesPosition,N,L,nlist)
                
                for i = 1:length(movedParticles)
                    movedP = movedParticles(i);
                    % recalculates pair distances after moving a particle    

                    xi = newParticlesPosition(1,movedP);
                    yi = newParticlesPosition(2,movedP);
                    newdist = dist;

                    % recalculate the relevent row elements in dist matrix

                    if movedP > 1
                        if ~isempty(nlist)
                            neiInd =...
                                nlist.neighborsindy(nlist.neighborsindx == movedP);
                            newdist(movedP,neiInd) =...
                            distPBC(xi,yi,...
                            newParticlesPosition(:,neiInd),...
                            L);

                        else

                            newdist(movedP,1:(movedP-1)) =...
                                distPBC(xi,yi,...
                                newParticlesPosition(:,1:(movedP-1)),L);
                        end
                    end

                    % recalculate the relevent column elements in dist matrix

                    if movedP < N
                        if ~isempty(nlist)
                            neiInd =...
                                nlist.neighborsindx(nlist.neighborsindy == movedP);
                            newdist(neiInd,movedP) =...
                            distPBC(xi,yi,...
                                newParticlesPosition(:,neiInd),L);

                        else
                            newdist((movedP + 1):N,movedP) =...
                            distPBC(xi,yi,...
                                newParticlesPosition(:,(movedP+1):N),L);
                        end
                    end
                end   
        end
    
        function dU = Uchange(movedParticle,dist,newDist,N,rCutoff,m)
        % calculates the change in energy after a particle has moved
                
                % calculate the old energy for the relevant particle pairs
                
                if movedParticle > 1
                    oldUrow = ...
                        pairU(dist(movedParticle,1:(movedParticle - 1)),rCutoff,m);
                else 
                    oldUrow = 0;
                end
                
                if movedParticle < N
                    oldUcol = ...
                        pairU(dist((movedParticle + 1):N,movedParticle),rCutoff,m);
                else 
                    oldUcol = 0;
                end

                oldU = oldUrow + oldUcol;
                
                % calculate the new energy for the relevant particle pairs
                
                if movedParticle > 1
                    newUrow = pairU(newDist...
                            (movedParticle,1:(movedParticle - 1)),rCutoff,m);
                else 
                    newUrow = 0;
                end
                
                if movedParticle < N
                    newUcol = pairU(newDist...
                        ((movedParticle + 1):N,movedParticle),rCutoff,m);
                else 
                    newUcol = 0;
                end
                
                newU = newUrow + newUcol;
                
                % clculate the change in energy
                
                dU = newU - oldU;
                   
        end
    
        function dV = Vchange(movedParticle,dist,newDist,N,rCutoff,rho,m)
        % calculates the change in the virial after a particle has moved
                
                % calculate the old virial for the relevant particle pairs
                
                if movedParticle > 1
                    oldVrow = ...
                        calcVirial(dist(movedParticle,1:(movedParticle - 1))...
                        ,rho,12,m,N,rCutoff);
                else 
                    oldVrow = 0;
                end
                
                if movedParticle < N
                    oldVcol = ...
                        calcVirial(dist((movedParticle + 1):N,movedParticle)...
                        ,rho,12,m,N,rCutoff);
                else 
                    oldVcol = 0;
                end

                oldV = oldVrow + oldVcol;
                
                % calculate the new virial for the relevant particle pairs
                
                if movedParticle > 1
                    newVrow = calcVirial(newDist...
                            (movedParticle,1:(movedParticle - 1))...
                            ,rho,12,m,N,rCutoff);
                else 
                    newVrow = 0;
                end
                
                if movedParticle < N
                    newVcol = calcVirial(newDist...
                        ((movedParticle + 1):N,movedParticle)...
                        ,rho,12,m,N,rCutoff);
                else 
                    newVcol = 0;
                end
                
                newV = newVrow + newVcol;
                
                % clculate the change in the virial
                
                dV = newV - oldV;
        end
        
    
   function dist = distPBC(x,y,allPossition,L)
            % calculate the PBC distances of the point (x,y) from all
            % particle possitions
            
                distx = abs(x - allPossition(1,:));
                disty = abs(y - allPossition(2,:));
                bigDistx = distx > (L/2);
                bigDisty = disty > (L/2);
                distx = distx - L.*bigDistx;
                disty = disty - L.*bigDisty;
                dist = sqrt(distx.^2 + disty.^2);
   end         
    function U = pairU(dist,rCutoff,m)

    % calculates the reduced energy according to the pair
    % potantial, only pair closer than rCutoff are regarded.

    % input: dist is a row vector of all pair distances
    % output: U is the total energy 

        dist_lt_rCutoff = dist(dist < rCutoff);
        % u is the energies of each pair
        u = 4*(((1./dist_lt_rCutoff).^12)-((1./dist_lt_rCutoff).^(-6))); 
        U = sum(u);
        
    
    end

    function virial = calcVirial(dists,rho,n,m,N,rCutoff)

    %% calculate the reduced virial of a system of 2D particles with pair potantial
    %       using the distances between particles 'allDist',
    %       calculate the virial.

    %        ==================================
    %       | virial = - (1/2)sum(r_ij (dU/dr)) |
    %        ==================================

    %       in reduced units, with a 4epsilon[(1/r)^n - (1/r)^m] potantial
    %       this is equivivalent to:
    %    
    %        ===========================================================
    %       | virial' = - (2*rho'/N)*sum(m(1/r'_ij)^m - n(1/r'_ij)^n) |
    %        ===========================================================

    %       where r_ij is the distance between particles i,j. the sum is on all
    %       particle pairs, each pair is counted once.
    %       for more information about the first equation: http://www2.msm.ctw.utwente.nl/sluding/TEACHING/APiE_Script_v2011.pdf
    %       eq. (3.15)       
    %       for more information about the second equation, read 'a note about
    %       reduced units' in the git repository

    %       input: 
    %       allDist - a matrix containing all pair distances for each step in a
    %                 monte carlo simulation, created with 'MonteCarlo2DLJ'
    %       rho - density in reduced units
    %       n - the pair potential power law in the distance 
    %           (pair potential in reduced units is u = 4*[(1/r)^n - (1/r)^m] 
    %           where r is the pair distance)
    %       m - the 'm' constant in the pair potantial
    %       T - reduced Temperature 

           dists = dists(and(dists < rCutoff,dists > 0));

           % V = - (2*rho'/N)*sum(m(1/r'_ij)^m - n(1/r'_ij)^n)
           vir = - (2*rho/N)*(sum(sum(6*dists.^(-6) - 12*dists.^(-12))));
           % make vir a row vector
           virial(1,:) = vir(1,1,:);

           
         
    end
end

