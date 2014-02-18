function [tout,yout,vout,rout,lout,stats] = NewmarkC3D(residual,yjac,constraints,tspan,y0,v0,l0)

%       [tout,yout,vout,rout,lout,stats] = NewmarkC3D(residual,yjac,constraints,tspan,y0,v0,l0)
%
%       Function solves constrained second order differential equations by Newmark/trapezoid-rule.
%       Function integrates the system of second order differential equations M(y)y'' = f_internal(t,y,v)
%       starting from time t0 to tfinal with initial conditions y0, v0, and fixed time step h.
%       Called function residual must return a column vector.
%
%       INPUT:
%       residual        function that calculates internal force vector
%       yjac            function that calculates mass, damping and stiffness matrices (jacobian matrices)
%       constraints     function that calculates constraint vector
%       tspan           time span vector, [t_start, t_finish, h]        (h=time step)
%       y0              initial position vector (use actual initial rotations)
%       v0              initial velocity vector
%       l0              initial Lagrangian multiplier vector (eg. zero vector)
%
%       OUTPUT:
%       tout            time vector
%       yout            displacement solution vector
%       vout            velocity solution vector
%       rout            acceleration solution vector
%       lout            Lagrangian multiplier vector
%       stats           statistic vector
%
%       NOTES           ID-table is required also for the rigid dynamic problems due to combability
%                       with flexible multibody problems. This code assumes six nodal degrees of freedom.
%                       For rigid elements there is one node per element. Additional masses do not have their
%                       own nodes
%
%       VERSION         17/02/2014
%
%      

% global variables

global FEM;
global Counter

% parameters
false = 0;              % logical variable
true = ~false;          % logical variable
mod_Newton = true;      % modified Newton-Raphson-iteration
MAX_ITERA=20000;          % max number of iteration per time step
epsilon=1e-3;           % residual norm convergence criterion
eta=1e-6;               % constraint convergence treshold

% Newmark parameters
t = tspan(1);           % initial time
tfinal = tspan(2);      % final time
h = tspan(3);           % time step
Alfa = 0.01;            % Damping parameter
Gamma = 1/2+Alfa;
Beta = 1/4*(Gamma+1/2)^2;
h2 = h*h;
hg = h*Gamma;
ophg = 1/hg;
ophg2 = ophg*ophg;
opbh2 = 1/(h*h*Beta);
gphbeta = Gamma/(h*Beta);

% statistics
nsteps = 0;             % stats
nfevals = 0;            % stats
nsolves = 0;            % stats
npds = 0;               % stats

% Initial acceleration +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

neq = length(y0);       % number of degrees of freedom
nconst=length(l0);      % number of constraints

% we collect the first row of output Note that the initial acceleration is zero
Counter = 1;
tout(Counter,1) = t;
yout(Counter,:) = y0.';
y1 = y0;

vout(Counter,:) = v0.';
rout(Counter,:) = zeros(neq,1)';
lout(Counter,:) = l0.';

%%% Residual invitation to the modified-Raphson method
feval(residual,t,y0,v0,v0*0,l0,0);
%%%

% for i=1:FEM.NumEle;
%     for j=1:FEM.ELE(i).NumNodes;
%         node=FEM.ELE(1,i).Node(j);
%         coo(node,:)=FEM.ELE(1,i).X((1+3*(j-1)):(3+3*(j-1)));
%     end
% end


% THE MAIN LOOP ***********************************************************************
done = false;
while ~done
    % LOOP FOR ADVANCING ONE STEP.
    iter = 0;
    
    % initial values for iteration
    
    

    y1 = y1 + h*vout(Counter,:)' + h2*(1/2-Beta)*rout(Counter,:)';
    
    %initial rotation is now set to zero
    
    for i=1:FEM.NumNodes
        for j=1:3
            if FEM.ID(i,j+3)>0
                y1(FEM.ID(i,j+3))=0;
            end % if
        end % for j
    end % for i
    
    v1 = vout(Counter,:)' + h*(1-Gamma)*rout(Counter,:)';
    a1 = zeros(neq,1);
    l1 = lout(Counter,:)';
    
    if mod_Newton % Modified Netwon-Raphson method
        % jacobian matrices
        [M,C,K,B] = feval(yjac,t+h,y1,v1,a1,l1);
        npds = npds + 1;                    % stats
        St = K + gphbeta*C + opbh2*M;
        A = [St, B'; B, zeros(nconst)];
        npds = npds + 1;                    % stats
    end
    
    nofailed = true;
    
    % internal iteration loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while true
        % residual vector
        r = feval(residual,t+h,y1,v1,a1,l1,0);
        % vector of costraints
        Psii=feval(constraints,t+h,y1,v1,a1);
        
        % jacobian matrices
        if ~mod_Newton  % Newton Raphson iteration. 
            [M,C,K,B] = feval(yjac,t+h,y1,v1,a1,l1);
            npds = npds + 1;                    % stats
            St = K + gphbeta*C + opbh2*M;
            A = [St,B';B,zeros(nconst)];
        end
        
        y = [r; -Psii];
        x = A\y;  % resulotion of the linearized system
        
        dq = x(1:neq,1);
        dl = x(neq+1:neq+nconst,1);
        
        iter = iter + 1;
        
        % rnorm = norm(r, inf);
        % psiinorm = norm(Psii, inf);
        rnorm = norm(r);
        psiinorm = norm(Psii);
        
        % update stage values
        l1 = l1 + dl;
        y1 = y1 + dq;
        v1 = v1 + gphbeta*dq;
        a1 = a1 + opbh2*dq;
        
        ymille=reshape(y1,6,FEM.NumNodes);
        ymille=ymille';
        ymille=ymille(:,1:3);
        
        displaced_beam=coo+ymille;
        
        % Newton-Raphson-convergence criterion
        if rnorm < epsilon && psiinorm < eta,
            break;
             elseif iter > MAX_ITERA,                % too many iterations ?
            done = true;
            disp('Newton-Raphson iteration failed !!!')
            nofailed = false;
            break;
        end %if
    end     % internal iteration loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    nfevals = nfevals + iter;           % stats
    nsolves = nsolves + iter;           % stats
    
    if nofailed,
        
        % Hiilihanko update mallin elementeille, joilla on kiertymävektori (massalla ei ole)
        y2 = y1;  %otetaan kiertymämuutokset ensin talteen
        r =  feval(residual,t+h,y1,v1,a1,l1,1);
        
        for i=1:FEM.NumNodes
            for j=1:3
                if FEM.ID(i,3+j)==0
                    Hiilihanko(j,1)=0;
                else
                    Hiilihanko(j,1)=y1(FEM.ID(i,3+j),1);
                end
            end
            
            
            R=MB_R(Hiilihanko);
            %FEM.Ele(i).Rref=FEM.Ele(i).R0*R;
            % sitten vielä nollataan kiertymäsuureet
            for j=1:3
                if FEM.ID(i,3+j)>0
                    y1(FEM.ID(i,3+j),1)=0;
                end
            end
        end
        
        nsteps = nsteps +1;             % stats
        Counter = Counter + 1;
        tout(Counter,1) = t+h;
        vout(Counter,:) = v1.';
        rout(Counter,:) = a1.';
        lout(Counter,:) = l1.';
        % sitten vielä siirtymäsuureet
        %yout(Counter,:) = y1.';
        for i=1:FEM.NumNodes
            jump=(i-1)*6;
            % translaatioliike
            for j=1:3
                yout(Counter,j+jump)=y2(FEM.ID(i,j));
            end
            % ja sitten rotaatio
            FEM.Kahveli(i).R=FEM.Kahveli(i).R*MB_R([y2(FEM.ID(i,4),1);y2(FEM.ID(i,5),1);y2(FEM.ID(i,6),1)]);
            Hiilihanko=MB_RTOY(FEM.Kahveli(i).R);
            yout(Counter,jump+4)=Hiilihanko(1,1);
            yout(Counter,jump+5)=Hiilihanko(2,1);
            yout(Counter,jump+6)=Hiilihanko(3,1);
        end
        
        if Counter==24
            t=t
        end
        t = t + h;                  % time incrementation
        
        if t >= tfinal,             % stop criterion
            done = true;
        end
        
        if mod(Counter,20),             % echo time to just see that something happens
            t
        end%
    end
end

% print cost statistics
fprintf('%g successful steps\n', nsteps);
fprintf('%g calls to functions\n', nfevals);
fprintf('%g calls to tangents \n', npds);
fprintf('%g solutions of linear systems\n', nsolves);

stats = [nsteps; nfevals; npds; nsolves];

end %function NewmarkC