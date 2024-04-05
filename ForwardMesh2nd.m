classdef ForwardMesh2nd < handle
%This is a second order forward-mesh class. NOTE: the Jacobian is still
%computed in first order basis.


    properties
        
        g   %the node points (an n x m array, each row containing one point)
        ng  %number of nodes
        gdim%dimension of the mesh
        H   %the connectivity matrix (each row contains node indices (row number of g) that are in a single element
        nH  %number of elements
        E   %Electrode surface elements, nel x 1 cell array, each cell n x gdim array, indices refer to rows in g
        nel %number of electrodes in the system
        EC  %Element connections; each row i contains all elements (rows of H) such that contain node i.
        itof%inverse-to-forward mesh transform matrix (1 if no inverse mesh is specified.)
        nginv%number of nodes in the inverse mesh
    end
    
    
    methods
        
        function obj = ForwardMesh2nd(g, H, E)
            
            obj.g = g;
            obj.H = H;
            obj.E = E;
            obj.gdim = size(g,2);
            
            obj.ng = size(g,1);
            obj.nH = size(H,1);
            obj.nel = length(E);
            obj.SetEC();
            obj.itof = 1;
            
        end
        
        function s1 = ItoF(self, sigma)
            s1 = self.itof*sigma;
        end
        
        function JI = JacobianFtoI(self, J)
            JI = J*self.itof;
        end
        
        function A = SigmadPhiidPhij(self, sigma)
           %This function returns matrix A containing in its element (i,j)
           %the integral of sigma*grad(phi_i) dot grad(phi_j).
            
            sigma = self.itof*sigma;

            k = 1;  
            if self.gdim == 2
                nphi = 6;
            elseif self.gdim == 3
                nphi = 10;
            end
            Arow = zeros(nphi*self.nH,nphi);
            Acol = zeros(nphi*self.nH,nphi);
            Aval = zeros(nphi*self.nH,nphi);

            % Gauss quadrature points and weights
            if self.gdim == 3
                a=0.58541020;b=0.13819660;
                ip=[b b b;a b b;b a b;b b a];
            elseif self.gdim == 2
                ip=[0.5 0; 0.5 0.5; 0 0.5];
            end

            L = cell(size(ip,1),1);
            for ii = 1:size(ip,1)
                if self.gdim == 3
                    L{ii}= [4*ip(ii,1)+4*ip(ii,2)+4*ip(ii,3)-3 4*ip(ii,1)-1 0 0 -8*ip(ii,1)+4-4*ip(ii,2)-4*ip(ii,3) 4*ip(ii,2) -4*ip(ii,2) -4*ip(ii,3) 4*ip(ii,3) 0;
                            4*ip(ii,1)+4*ip(ii,2)+4*ip(ii,3)-3 0 4*ip(ii,2)-1 0 -4*ip(ii,1) 4*ip(ii,1) -8*ip(ii,2)+4-4*ip(ii,1)-4*ip(ii,3) -4*ip(ii,3) 0 4*ip(ii,3);
                            4*ip(ii,1)+4*ip(ii,2)+4*ip(ii,3)-3 0 0 4*ip(ii,3)-1 -4*ip(ii,1) 0 -4*ip(ii,2) -8*ip(ii,3)+4-4*ip(ii,1)-4*ip(ii,2) 4*ip(ii,1) 4*ip(ii,2)];
                elseif self.gdim == 2
                    L{ii} = error;
                end
            end
            for ii=1:self.nH
              % Go through all triangles/tetrahedra
              ind = self.H(ii,:);
              gg = self.g(ind,:);
              ss = sigma(ind(1:self.gdim+1));
              int = self.intLinSigma(gg,ss,ip,L);
              Arow(k:k+nphi-1,:) = repmat(ind', 1, nphi);
              Acol(k:k+nphi-1,:) = repmat(ind, nphi, 1);
              Aval(k:k+nphi-1,:) = int; 

              k = k + nphi;
            end  
            
            A = sparse(Arow,Acol,Aval,self.ng+self.nel-1,self.ng+self.nel-1);

            
        end
        
        function [S, M, B] = EITElectrodeTerms(self)
            %This function returns matrices S and M, and vector B.
            %S is a cell array containing integrals phi_i*phi_j along
            %each electrode surfaces (one electrode per cell)
            %M contains integral phi_i on surface of electrode j
            %B contains the electrode areas in its elements.

            B = zeros(self.nel,1);

            M = zeros(self.ng, self.nel);

            S = cell(self.nel,1);

            if self.gdim == 2
                nphi = 3;
            elseif self.gdim == 3
                nphi = 6;
            end
            % Loop through electrodes
            for ii=1:self.nel
              spos = 1;
              faces = self.E{ii};
              len = self.gdim*size(self.E{ii},1);
              intS = zeros(len,nphi);
              rowS = zeros(len,nphi);
              colS = zeros(len,nphi);
              % Loop through faces on electrode ii
              for jj = 1:size(faces,1)
                ind = faces(jj,:); % face nodes

                nnodes = length(ind);
                gg = self.g(ind,:);
                rcidx = repmat(ind, nnodes, 1);

                if self.gdim == 3
                    bb1 = self.triang1(gg);
                    bb2 = self.triang2(gg);
                elseif self.gdim == 2
                    error('2nd order 2D mesh not yet supported')
                    bb1 = self.phii1D(gg);
                    bb2 = self.phiiphij1D(gg);
                end

                intS(spos:spos+nnodes-1,:) = bb2;
                rowS(spos:spos+nnodes-1,:) = rcidx.'; 
                colS(spos:spos+nnodes-1,:) = rcidx;   

                M(ind,ii) = M(ind,ii) + bb1;

                if self.gdim == 3
                    B(ii,:)  = B(ii,:) + self.elektro(gg);
                elseif self.gdim == 2
                    B(ii,:)  = B(ii,:) + self.ElectrodeArea1D(gg);
                end

                spos = spos + nnodes;
              end
              S{ii} = sparse(rowS,colS,intS,self.ng,self.ng);
            end
            
        end
        
        function [Ai, Av] = GradientMatrix(self)
            %Computes all the basis function gradients of the 1st order
            %mesh to be used in calculating the Jacobian.
            ar = zeros(self.ng*12,self.gdim+1);
            ac = zeros(self.ng*12,self.gdim+1);
            av = zeros(self.ng*12,self.gdim+1);

            Ai = cell(self.ng,1);
            Av = cell(self.ng,1);

            % compute gradients
            for jj=1:self.ng

              El = self.EC(jj,self.EC(jj,:)>0);

              rid = 1;
              for ii=1:length(El)
                ind = self.H(El(ii),1:self.gdim+1); % Indices of the element
                gg=self.g(ind,:);

                idc = repmat(ind,self.gdim+1,1);
                idr = idc';
                
                L=[-ones(self.gdim,1) eye(self.gdim)];
                Jt=L*gg;
                dJt=abs(det(Jt)); % Tetrahedra volume
                G=Jt\L; % Gradients of each basis function
                GdJt=G'*G*dJt;

                if self.gdim == 3
                    int=1/24*GdJt;
                elseif self.gdim == 2
                    int=1/6*GdJt;
                end

                % temporary storage
                ar(rid:rid+self.gdim,:) = idr;
                ac(rid:rid+self.gdim,:) = idc;
                av(rid:rid+self.gdim,:) = int;
                rid = rid + 1 + self.gdim;      
              end     
              I = 1:(rid-1);

              S = sparse(ar(I,:),ac(I,:),av(I,:),self.ng,self.ng);
              [I,~] = find(S);

              L = unique(I);    
              Ai{jj} = L;             % symmetric assumption => no columns needed
              Av{jj} = full(S(L,L));
            end
        end
        
        function SetInverseMesh(self, imesh)
            
            if size(imesh.g,2) == 2%Check the dimension of the inverse mesh
                self.itof = self.interpolatematrix2d(imesh.H, imesh.g, self.g(:,1:2));
            elseif size(imesh.g,2) == 3
                self.itof = interpolatematrix3d(imesh.H, imesh.g, self.g);
            else
                error(['Unfit second dimension of ginv: size(ginv,2) = ' num2str(size(ginv,2))]);
            end
            self.nginv = size(imesh.g,1);
            
        end
        
        function SetEC(self)
            %Sets self.EC, i.e. an array, for which each row i contains all
            %such element indices (i.e. rows of H) that the element contains the node i (i.e. row i of g).
            %currently the only use for self.EC is by self.GradientMatrix()
            counter = zeros(self.ng, 1);
            self.EC = zeros(self.ng, 4);
            for iH = 1:self.nH
                for in = 1:size(self.H,2)
                    id = self.H(iH, in);
                    counter(id) = counter(id) + 1;
                    self.EC(id, counter(id)) = iH;
                end
            end
            
        end
        
        function int=intLinSigma(self,g,s,ip,L)
            S = [1-sum(ip,2) ip]*s;
            int = 0;
            for ii = 1:length(L)
                Jt = L{ii}*g;
                iJt = inv(Jt);
                dJt = abs(det(Jt));
                G = iJt*L{ii};
                GdJt = G'*G*dJt/factorial(self.gdim+1);
                int = int + S(ii)*GdJt;
            end

        end

    end

    methods(Static)

        function int=triang1(g)

            w=[1/40;1/15;1/40;1/15;1/40;1/15;9/40];
            ip=[0 0; 1/2 0;1 0;1/2 1/2;0 1;0 1/2;1/3 1/3];
            
            x=g(:,1);
            y=g(:,2);
            z=g(:,3);
            
            int=0;
            for ii=1:length(ip)
            S=[2*ip(ii,1)^2+2*ip(ii,2)^2+4*ip(ii,1)*ip(ii,2)-3*ip(ii,1)-3*ip(ii,2)+1; ...
               2*ip(ii,1)^2-ip(ii,1); ...
               2*ip(ii,2)^2-ip(ii,2); ...
               -4*ip(ii,1)^2-4*ip(ii,1)*ip(ii,2)+4*ip(ii,1); ...
               4*ip(ii,1)*ip(ii,2); ...
               -4*ip(ii,2)^2-4*ip(ii,1)*ip(ii,2)+4*ip(ii,2)];
            
            k=ip(ii,1);
            n=ip(ii,2);
            c=0;
            
            aa=2*(2*x(1)+2*x(2)-4*x(4))*k+(-3*x(1)-x(2)+4*x(4))+(4*x(1)-4*x(4)+4*x(5)-4*x(6))*n;
            bb=2*(2*y(1)+2*y(2)-4*y(4))*k+(-3*y(1)-y(2)+4*y(4))+(4*y(1)-4*y(4)+4*y(5)-4*y(6))*n;
            cc=2*(2*z(1)+2*z(2)-4*z(4))*k+(-3*z(1)-z(2)+4*z(4))+(4*z(1)-4*z(4)+4*z(5)-4*z(6))*n;
            AA=2*(2*x(1)+2*x(3)-4*x(6))*n+(-3*x(1)-x(3)+4*x(6))+(4*x(1)-4*x(4)+4*x(5)-4*x(6))*k;
            BB=2*(2*y(1)+2*y(3)-4*y(6))*n+(-3*y(1)-y(3)+4*y(6))+(4*y(1)-4*y(4)+4*y(5)-4*y(6))*k;
            CC=2*(2*z(1)+2*z(3)-4*z(6))*n+(-3*z(1)-z(3)+4*z(6))+(4*z(1)-4*z(4)+4*z(5)-4*z(6))*k;
            
            ll=sqrt((bb*CC-BB*cc)^2+(cc*AA-CC*aa)^2+(aa*BB-AA*bb)^2);
            
            int=int+w(ii)*ll*S;
             
            end
        end
        
        function int=triang2(g)
            
            w=[1/40;1/15;1/40;1/15;1/40;1/15;9/40];
            ip=[0 0; 1/2 0;1 0;1/2 1/2;0 1;0 1/2;1/3 1/3];
            
            x=g(:,1);
            y=g(:,2);
            z=g(:,3);
            
            int=0;
            for ii=1:length(ip)
            S=[2*ip(ii,1)^2+2*ip(ii,2)^2+4*ip(ii,1)*ip(ii,2)-3*ip(ii,1)-3*ip(ii,2)+1; ...
               2*ip(ii,1)^2-ip(ii,1); ...
               2*ip(ii,2)^2-ip(ii,2); ...
               -4*ip(ii,1)^2-4*ip(ii,1)*ip(ii,2)+4*ip(ii,1); ...
               4*ip(ii,1)*ip(ii,2); ...
               -4*ip(ii,2)^2-4*ip(ii,1)*ip(ii,2)+4*ip(ii,2)];
            
            k=ip(ii,1);
            n=ip(ii,2);
            c=0;
            
            aa=2*(2*x(1)+2*x(2)-4*x(4))*k+(-3*x(1)-x(2)+4*x(4))+(4*x(1)-4*x(4)+4*x(5)-4*x(6))*n;
            bb=2*(2*y(1)+2*y(2)-4*y(4))*k+(-3*y(1)-y(2)+4*y(4))+(4*y(1)-4*y(4)+4*y(5)-4*y(6))*n;
            cc=2*(2*z(1)+2*z(2)-4*z(4))*k+(-3*z(1)-z(2)+4*z(4))+(4*z(1)-4*z(4)+4*z(5)-4*z(6))*n;
            AA=2*(2*x(1)+2*x(3)-4*x(6))*n+(-3*x(1)-x(3)+4*x(6))+(4*x(1)-4*x(4)+4*x(5)-4*x(6))*k;
            BB=2*(2*y(1)+2*y(3)-4*y(6))*n+(-3*y(1)-y(3)+4*y(6))+(4*y(1)-4*y(4)+4*y(5)-4*y(6))*k;
            CC=2*(2*z(1)+2*z(3)-4*z(6))*n+(-3*z(1)-z(3)+4*z(6))+(4*z(1)-4*z(4)+4*z(5)-4*z(6))*k;
            
            ll=sqrt((bb*CC-BB*cc)^2+(cc*AA-CC*aa)^2+(aa*BB-AA*bb)^2);
            
            int=int+w(ii)*ll*S*S';
             
            end


        end
        
        function int=phii1D(g)
            
            a = 1/2;
            l=sqrt((g(1,1)-g(2,1))^2 + (g(1,2)-g(2,2))^2);
            int = l*[a;a];

        end

        function int=phiiphij1D(g)
            
            a = 1/6;
            b = 1/3;
            l=sqrt((g(1,1)-g(2,1))^2 + (g(1,2)-g(2,2))^2);
            int=l*[b a;a b]; 

        end

        function int=elektro(g)
            w=[1/40;1/15;1/40;1/15;1/40;1/15;9/40];
            ip=[0 0; 1/2 0;1 0;1/2 1/2;0 1;0 1/2;1/3 1/3];
            
            x=g(:,1);
            y=g(:,2);
            z=g(:,3);
            
            int=0;
            for ii=1:length(ip)
            
            k=ip(ii,1);
            n=ip(ii,2);
            c=0;
            aa=2*(2*x(1)+2*x(2)-4*x(4))*k+(-3*x(1)-x(2)+4*x(4))+(4*x(1)-4*x(4)+4*x(5)-4*x(6))*n;
            bb=2*(2*y(1)+2*y(2)-4*y(4))*k+(-3*y(1)-y(2)+4*y(4))+(4*y(1)-4*y(4)+4*y(5)-4*y(6))*n;
            cc=2*(2*z(1)+2*z(2)-4*z(4))*k+(-3*z(1)-z(2)+4*z(4))+(4*z(1)-4*z(4)+4*z(5)-4*z(6))*n;
            AA=2*(2*x(1)+2*x(3)-4*x(6))*n+(-3*x(1)-x(3)+4*x(6))+(4*x(1)-4*x(4)+4*x(5)-4*x(6))*k;
            BB=2*(2*y(1)+2*y(3)-4*y(6))*n+(-3*y(1)-y(3)+4*y(6))+(4*y(1)-4*y(4)+4*y(5)-4*y(6))*k;
            CC=2*(2*z(1)+2*z(3)-4*z(6))*n+(-3*z(1)-z(3)+4*z(6))+(4*z(1)-4*z(4)+4*z(5)-4*z(6))*k;
            
            ll=sqrt((bb*CC-BB*cc)^2+(cc*AA-CC*aa)^2+(aa*BB-AA*bb)^2);
            
            int=int+w(ii)*ll;
             
            end

    
        end

        function int=ElectrodeArea1D(g)

            int = sqrt((g(1,1)-g(2,1))^2 + (g(1,2)-g(2,2))^2);
        
        end

        function M = interpolatematrix2d(H,g,p)
            %Use this: https://mathworld.wolfram.com/TriangleInterior.html
            v0 = g(H(:,1),:);
            v1 = g(H(:,2),:) - v0;
            v2 = g(H(:,3),:) - v0;
            a = (p(:,1)*v2(:,2)' - p(:,2)*v2(:,1)' - repmat((v0(:,1).*v2(:,2)-v0(:,2).*v2(:,1))',size(p,1),1))./repmat((v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1))',size(p,1),1);
            b = -(p(:,1)*v1(:,2)' - p(:,2)*v1(:,1)' - repmat((v0(:,1).*v1(:,2)-v0(:,2).*v1(:,1))',size(p,1),1))./repmat((v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1))',size(p,1),1);
            insides = a>0 & b>0 & a+b<1-1e-12;
            T = zeros(size(p,1), 1);
            for ip = 1:size(p,1)
                temp = find(insides(ip,:));
                if length(temp) > 1
                    error('a point is inside several triangles!');
                elseif isempty(temp)
                    T(ip) = nan;
                else
                    T(ip) = temp;
                end
            end
            nT = length(T);
            
            p = p';
            L = [-1 -1;1 0;0 1];
            f = zeros(nT,1);
            
            M = sparse(size(p,2),size(g,1));
            
            for ii=1:nT
              iH = T(ii);
              
              if ~isnan(iH)
                id = H(iH,:);
                gg = g(id,:)';
              
                X = (gg*L)\(p(:,ii)-gg(:,1));
            
                M(ii,id) = M(ii,id) + [1-X(1)-X(2), X(1), X(2)];
              else
                dist_to_nearest = 1e9;
                nearest = 0;
                for j=1:size(g,1)
                    dist = sqrt((p(1,ii)-g(j,1)).^2+(p(2,ii)-g(j,2)).^2);
                    if(dist<dist_to_nearest)
                        dist_to_nearest = dist;
                        nearest = j;
                    end
                end
                M(ii,nearest) = M(ii,nearest) + 1;
              end  
            end
        end
        
        function M = interpolatematrix3d(H,g,p)
            %Use this: http://steve.hollasch.net/cgindex/geometry/ptintet.html
            
            p1 = g(H(:,1),:);
            p2 = g(H(:,2),:);
            p3 = g(H(:,3),:);
            p4 = g(H(:,4),:);
            
            D0 = p2(:,1).*p3(:,3).*p4(:,2) + p2(:,2).*p3(:,1).*p4(:,3) + p2(:,3).*p3(:,2).*p4(:,1) ...
                -p2(:,1).*p3(:,2).*p4(:,3) - p2(:,2).*p3(:,3).*p4(:,1) - p2(:,3).*p3(:,1).*p4(:,2) ...
                -p1(:,1).*p3(:,3).*p4(:,2) - p1(:,2).*p3(:,1).*p4(:,3) - p1(:,3).*p3(:,2).*p4(:,1) ...
                +p1(:,1).*p3(:,2).*p4(:,3) + p1(:,2).*p3(:,3).*p4(:,1) + p1(:,3).*p3(:,1).*p4(:,2) ...
                +p1(:,1).*p2(:,3).*p4(:,2) + p1(:,2).*p2(:,1).*p4(:,3) + p1(:,3).*p2(:,2).*p4(:,1) ...
                -p1(:,1).*p2(:,2).*p4(:,3) - p1(:,2).*p2(:,3).*p4(:,1) - p1(:,3).*p2(:,1).*p4(:,2) ...
                -p1(:,1).*p2(:,3).*p3(:,2) - p1(:,2).*p2(:,1).*p3(:,3) - p1(:,3).*p2(:,2).*p3(:,1) ...
                +p1(:,1).*p2(:,2).*p3(:,3) + p1(:,2).*p2(:,3).*p3(:,1) + p1(:,3).*p2(:,1).*p3(:,2);
            
            D1 = p2(:,1).*p3(:,3).*p4(:,2) + p2(:,2).*p3(:,1).*p4(:,3) + p2(:,3).*p3(:,2).*p4(:,1) ...
                -p2(:,1).*p3(:,2).*p4(:,3) - p2(:,2).*p3(:,3).*p4(:,1) - p2(:,3).*p3(:,1).*p4(:,2) ...
                -p3(:,3).*p4(:,2)*p(:,1)' - p3(:,1).*p4(:,3)*p(:,2)' - p3(:,2).*p4(:,1)*p(:,3)' ...
                +p3(:,2).*p4(:,3)*p(:,1)' + p3(:,3).*p4(:,1)*p(:,2)' + p3(:,1).*p4(:,2)*p(:,3)' ...
                +p2(:,3).*p4(:,2)*p(:,1)' + p2(:,1).*p4(:,3)*p(:,2)' + p2(:,2).*p4(:,1)*p(:,3)' ...
                -p2(:,2).*p4(:,3)*p(:,1)' - p2(:,3).*p4(:,1)*p(:,2)' - p2(:,1).*p4(:,2)*p(:,3)' ...
                -p2(:,3).*p3(:,2)*p(:,1)' - p2(:,1).*p3(:,3)*p(:,2)' - p2(:,2).*p3(:,1)*p(:,3)' ...
                +p2(:,2).*p3(:,3)*p(:,1)' + p2(:,3).*p3(:,1)*p(:,2)' + p2(:,1).*p3(:,2)*p(:,3)';
            
            D2 = p3(:,3).*p4(:,2)*p(:,1)' + p3(:,1).*p4(:,3)*p(:,2)' + p3(:,2).*p4(:,1)*p(:,3)' ...
                -p3(:,2).*p4(:,3)*p(:,1)' - p3(:,3).*p4(:,1)*p(:,2)' - p3(:,1).*p4(:,2)*p(:,3)' ...
                -p1(:,1).*p3(:,3).*p4(:,2) - p1(:,2).*p3(:,1).*p4(:,3) - p1(:,3).*p3(:,2).*p4(:,1) ...
                +p1(:,1).*p3(:,2).*p4(:,3) + p1(:,2).*p3(:,3).*p4(:,1) + p1(:,3).*p3(:,1).*p4(:,2) ...
                +p1(:,1).*p4(:,2)*p(:,3)' + p1(:,2).*p4(:,3)*p(:,1)' + p1(:,3).*p4(:,1)*p(:,2)' ...
                -p1(:,1).*p4(:,3)*p(:,2)' - p1(:,2).*p4(:,1)*p(:,3)' - p1(:,3).*p4(:,2)*p(:,1)' ...
                -p1(:,1).*p3(:,2)*p(:,3)' - p1(:,2).*p3(:,3)*p(:,1)' - p1(:,3).*p3(:,1)*p(:,2)' ...
                +p1(:,1).*p3(:,3)*p(:,2)' + p1(:,2).*p3(:,1)*p(:,3)' + p1(:,3).*p3(:,2)*p(:,1)';
            
            D3 = p2(:,1).*p4(:,2)*p(:,3)' + p2(:,2).*p4(:,3)*p(:,1)' + p2(:,3).*p4(:,1)*p(:,2)' ...
                -p2(:,1).*p4(:,3)*p(:,2)' - p2(:,2).*p4(:,1)*p(:,3)' - p2(:,3).*p4(:,2)*p(:,1)' ...
                -p1(:,1).*p4(:,2)*p(:,3)' - p1(:,2).*p4(:,3)*p(:,1)' - p1(:,3).*p4(:,1)*p(:,2)' ...
                +p1(:,1).*p4(:,3)*p(:,2)' + p1(:,2).*p4(:,1)*p(:,3)' + p1(:,3).*p4(:,2)*p(:,1)' ...
                +p1(:,1).*p2(:,3).*p4(:,2) + p1(:,2).*p2(:,1).*p4(:,3) + p1(:,3).*p2(:,2).*p4(:,1) ...
                -p1(:,1).*p2(:,2).*p4(:,3) - p1(:,2).*p2(:,3).*p4(:,1) - p1(:,3).*p2(:,1).*p4(:,2) ...
                -p1(:,1).*p2(:,3)*p(:,2)' - p1(:,2).*p2(:,1)*p(:,3)' - p1(:,3).*p2(:,2)*p(:,1)' ...
                +p1(:,1).*p2(:,2)*p(:,3)' + p1(:,2).*p2(:,3)*p(:,1)' + p1(:,3).*p2(:,1)*p(:,2)';
            
            D4 = p2(:,1).*p3(:,3)*p(:,2)' + p2(:,2).*p3(:,1)*p(:,3)' + p2(:,3).*p3(:,2)*p(:,1)' ...
                -p2(:,1).*p3(:,2)*p(:,3)' - p2(:,2).*p3(:,3)*p(:,1)' - p2(:,3).*p3(:,1)*p(:,2)' ...
                -p1(:,1).*p3(:,3)*p(:,2)' - p1(:,2).*p3(:,1)*p(:,3)' - p1(:,3).*p3(:,2)*p(:,1)' ...
                +p1(:,1).*p3(:,2)*p(:,3)' + p1(:,2).*p3(:,3)*p(:,1)' + p1(:,3).*p3(:,1)*p(:,2)' ...
                +p1(:,1).*p2(:,3)*p(:,2)' + p1(:,2).*p2(:,1)*p(:,3)' + p1(:,3).*p2(:,2)*p(:,1)' ...
                -p1(:,1).*p2(:,2)*p(:,3)' - p1(:,2).*p2(:,3)*p(:,1)' - p1(:,3).*p2(:,1)*p(:,2)' ...
                -p1(:,1).*p2(:,3).*p3(:,2) - p1(:,2).*p2(:,1).*p3(:,3) - p1(:,3).*p2(:,2).*p3(:,1) ...
                +p1(:,1).*p2(:,2).*p3(:,3) + p1(:,2).*p2(:,3).*p3(:,1) + p1(:,3).*p2(:,1).*p3(:,2);
            
            
            sg1 = D0.*D1;
            sg2 = D0.*D2;
            sg3 = D0.*D3;
            sg4 = D0.*D4;
            
            tr = sg1>1e-12 & sg2>1e-12 & sg3>1e-12 & sg4>1e-12;
            
            T = zeros(size(p,1), 1);
            for ip = 1:size(p,1)
                s = sum(tr(:,ip));
                if s > 1
                    error('a point is inside several tetrahedra!');
                elseif s == 0
                    T(ip) = nan;
                else
                    T(ip) = find(tr(:,ip));
                end
            end
            
            p = p';
            L = [-1 -1 -1;1 0 0;0 1 0; 0 0 1];
            
            M = sparse(size(p,2),size(g,1));
            
            for ii=1:size(p,2)
              iH = T(ii);
              
              if ~isnan(iH)
                id = H(iH,:);
                gg = g(id,:)';
              
                X = (gg*L)\(p(:,ii)-gg(:,1));
            
                M(ii,id) = M(ii,id) + [1-X(1)-X(2)-X(3), X(1), X(2), X(3)];
              else
                dist_to_nearest = 1e9;
                nearest = 0;
                for j=1:size(g,1)
                    dist = sqrt((p(1,ii)-g(j,1)).^2+(p(2,ii)-g(j,2)).^2+(p(3,ii)-g(j,3)).^2);
                    if(dist<dist_to_nearest)
                        dist_to_nearest = dist;
                        nearest = j;
                    end
                end
                M(ii,nearest) = M(ii,nearest) + 1;
              end  
            end


    end
    
    
    
    end
    
    
    
    
    
    
    
end