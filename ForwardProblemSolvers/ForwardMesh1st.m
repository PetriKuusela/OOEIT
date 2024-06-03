classdef ForwardMesh1st < handle
%This is a first order forward-mesh class that can be used by the forward 
%problem solvers, e.g. EITFEM. It contains functions to compute
%all the mesh related things for the FEM, most notably different integrals
%of the basis functions, which are computed in functions SigmadPhiidPhij,
%EITElectrodeTerms and GradientMatrix.
%
%For implementing your own mesh class, in addition to the abovementioned
%functions, also functions ItoF, JacobianFtoI and SetInverseMesh are
%required, as these are called by the forward problem solver.


    properties
        
        g   %the node points (an n x m array, n = number of nodes, m = dimension of the mesh, each row contains coordinates of one point)
        ng  %number of nodes
        gdim%dimension of the mesh
        H   %the connectivity matrix (each row contains node indices (row number of g) that are in a single element)
        nH  %number of elements
        E   %Electrode surface elements, nel x 1 cell array, each cell n x gdim array (n = number of surface elements in that electrode), indices refer to rows in g
        nel %number of electrodes in the system
        EC  %Element connections; each row i contains all elements (rows of H) that contain node i.
        itof%inverse-to-forward mesh transform matrix (1 if no inverse mesh is specified.)
        nginv%number of nodes in the inverse mesh (if no inverse mesh is specified, it is the same as the forward mesh)
    end
    
    
    methods
        
        function obj = ForwardMesh1st(g, H, E)
            %Class constructor. For descriptions of the properties, see the
            %properties-section.

            obj.g = g;
            obj.H = H;
            obj.E = E;

            obj.gdim = size(g,2);
            obj.ng = size(g,1);
            obj.nH = size(H,1);
            obj.nel = length(E);
            obj.SetEC();
            obj.itof = 1;
            obj.nginv = obj.ng;
            
        end
        
        function s1 = ItoF(self, sigma)
            %This function transform sigma from inverse mesh basis to the
            %forward mesh basis. Inverse mesh is a separate mesh that can
            %be defined for the conductivity to reduce the number of
            %parameters to be estimated.
            s1 = self.itof*sigma;
        end
        
        function JI = JacobianFtoI(self, J)
            %This function implement the transformation of Jacobian from
            %the forward basis to the inverse mesh basis. The
            %reason for this multiplication comes from using the chain rule
            %of derivatives.
            JI = J*self.itof;
        end
        
        function A = SigmadPhiidPhij(self, sigma)
           %This function returns matrix A containing in its element (i,j)
           %the integral of sigma*grad(phi_i) dot grad(phi_j). That is also the 
           %conductivity dependent part of the upper left block of the
           %FEM matrix.
           %The input "sigma", is assumed to be in the inverse mesh basis.
            
            %Transform to forward mesh basis
            sigma = self.ItoF(sigma);

            %initialize variables where to collect the integral values. The
            %values are stored in row/column/value format and stored in a
            %sparse matrix in the end
            k = 1;  
            Arow = zeros((self.gdim+1)*self.nH,(self.gdim+1));
            Acol = zeros((self.gdim+1)*self.nH,(self.gdim+1));
            Aval = zeros((self.gdim+1)*self.nH,(self.gdim+1));

            % Gauss quadrature points for integration
            if self.gdim == 3
                a=0.58541020;b=0.13819660;
                ip=[b b b;a b b;b a b;b b a];
            elseif self.gdim == 2
                ip=[0.5 0; 0.5 0.5; 0 0.5];
            end

            % difference matrix of the (linear) basis functions
            L = [-ones(self.gdim,1) eye(self.gdim)];
            for ii=1:self.nH
              % Go through all triangles/tetrahedra
              ind = self.H(ii,:); %The indices of the nodes in element ii
              gg = self.g(ind,:); %The coordinates of the nodes in ind
              ss = sigma(ind);    %The conductivities at the nodes in ind

              %Integral of sigma*grad(phi_i) dot grad(phi_j) inside element ii:
              int = self.intLinSigma(gg,ss,ip,L); 
              %int is now a 3-by-3 or 4-by-4 matrix, with element i,j
              %containing the said integral with the basisfunctions phi_i
              %and phi_j being the basis functions of nodes ind(i) and
              %ind(j).
              Arow(k:k+self.gdim,:) = repmat(ind', 1, self.gdim+1);%the row indices of the int-values in the final sparse matrix to be constructed
              Acol(k:k+self.gdim,:) = repmat(ind, self.gdim+1, 1);%the column indices of the int-values in the final sparse matrix to be constructed
              Aval(k:k+self.gdim,:) = int; %the integral values

              k = k + self.gdim + 1;
            end  
            
            %The output:
            A = sparse(Arow,Acol,Aval,self.ng+self.nel-1,self.ng+self.nel-1);

            
        end
        
        function [S, M, B] = EITElectrodeTerms(self)
            %This function returns matrices S and M, and vector B.
            %S is a cell array containing integrals phi_i*phi_j along
            %each electrode surfaces (one electrode per cell)
            %M contains integral phi_i on surface of electrode j
            %B contains the electrode areas in its elements.

            B = zeros(self.nel,1);
            %Shouldn't this be sparse?
            M = zeros(self.ng, self.nel);
            S = cell(self.nel,1);

            % Loop through electrodes
            for ii=1:self.nel
              spos = 1;%index for storing values for matrix in cell array S
              faces = self.E{ii};%the boundary elements of electrode ii
              len = self.gdim*size(self.E{ii},1);%The number of rows in the matrices for collecting values for S-matrix
              intS = zeros(len,self.gdim);
              rowS = zeros(len,self.gdim);
              colS = zeros(len,self.gdim); %These three matrices contain the information for making the sparse matrix
              % Loop through faces on electrode ii
              for jj = 1:size(faces,1)
                ind = faces(jj,:); % face nodes
                gg = self.g(ind,:);% the coordinates of the nodes
                rcidx = repmat(ind, self.gdim, 1);

                %The following compute integrals of phi_i * phi_j (in bb2) 
                %and phi_i (in bb1) inside the boundary element jj
                if self.gdim == 3
                    [bb1, bb2] = self.phiiphij2D(gg);
                elseif self.gdim == 2
                    [bb1, bb2] = self.phiiphij1D(gg);
                end

                %Store the values
                intS(spos:spos+self.gdim-1,:) = bb2;
                rowS(spos:spos+self.gdim-1,:) = rcidx.'; 
                colS(spos:spos+self.gdim-1,:) = rcidx;   
                M(ind,ii) = M(ind,ii) + bb1;

                %compute the element area
                if self.gdim == 3
                    B(ii,:)  = B(ii,:) + self.ElementArea2D(gg);
                elseif self.gdim == 2
                    B(ii,:)  = B(ii,:) + self.ElementArea1D(gg);
                end

                spos = spos + self.gdim;
              end
              %All the elements of electrode ii have been gone through, so
              %next create the sparse matrix S{ii}
              S{ii} = sparse(rowS,colS,intS,self.ng,self.ng);
            end
            
        end
        
        function [Ai, Av] = GradientMatrix(self)
            %Computes all the basis function gradients of the
            %mesh. These are used in calculating the Jacobian.
            %
            %Ai and Av are cell arrays containing at cell jj the
            %information for the matrix obtained by differentiating the 
            %FEM matrix w.r.t the conductivity of node jj. Ai{jj} contains
            %the row/column indices (they're the same, as dA/dsigma is
            %symmetric) of non-zero elements of the derivative. Av{jj}
            %contains the corresponding values of integrals of grad(phi_i) dot
            %grad(phi_j), i.e. the entries of the derivative matrix.

            Ai = cell(self.ng,1);
            Av = cell(self.ng,1);

            for jj=1:self.ng
                %compute the derivative of the FEM matrix w.r.t
                %conductivity of node jj.

                %get the elements that contain the support of basis 
                %function jj, i.e. have the non-zero gradient dot products
                %in the derivative of the FEM matrix:
              El = self.EC(jj,self.EC(jj,:)>0);
              %due to how this info is stored (in a matrix), there may be
              %zeros on some (many) rows. They are to be omitted.

              %Initialize the matrices for collecting the indices and values
              %for creating a sparse matrix.
              ar = zeros(length(El)*(self.gdim+1),self.gdim+1);
              ac = zeros(length(El)*(self.gdim+1),self.gdim+1);
              av = zeros(length(El)*(self.gdim+1),self.gdim+1);

              %index for collecting values to the matrices
              rid = 1;
              for ii=1:length(El)%Go through the relevant elements
                  
                ind = self.H(El(ii),:); % Indices of the element
                gg=self.g(ind,:);%the node coordinates

                idc = repmat(ind,self.gdim+1,1);%The row and column indices of the values in the derivative matrix
                idr = idc';
                
                % difference matrix of the (linear) basis functions
                L=[-ones(self.gdim,1) eye(self.gdim)];
                Jt=L*gg;
                dJt=abs(det(Jt)); % Triangle/tetrahedra volume
                G=Jt\L; % Gradients of each basis function
                GdJt=G'*G*dJt;
                if self.gdim == 3
                    int=1/24*GdJt;
                elseif self.gdim == 2
                    int=1/6*GdJt;
                end
                %now int contains integrals of grad(phi_i) dot grad(phi_j)
                %for i and j in nodes in ind

                % temporary storage
                ar(rid:rid+self.gdim,:) = idr;
                ac(rid:rid+self.gdim,:) = idc;
                av(rid:rid+self.gdim,:) = int;
                rid = rid + 1 + self.gdim;      
              end     

              %Create the sparse derivative matrix of the FEM matrix
              S = sparse(ar,ac,av,self.ng,self.ng);
              %find non-zero rows/columns
              [I,~] = find(S);

              L = unique(I);  %remove duplicates  
              % symmetric assumption => no columns needed
              Ai{jj} = L;     %store the indices
              Av{jj} = full(S(L,L));%store the values located at the indices
            end
        end
        
        function SetInverseMesh(self, imesh)
            %A separate mesh can be used for representing the conductivity.
            %Here we call that mesh an inverse mesh. In this mesh object
            %(which stores the mesh for FEM computations of the forward
            %problem), we do not store the inverse mesh. Instead, we only
            %store the transformation matrix, that allows us to interpolate
            %values from inverse mesh to the forward mesh.
            %
            %This function computes and stores that matrix, based on given
            %inverse mesh "imesh", which should have two fields: g and H
            %(defined similarly as the g and H of this class)
            
            %The interploation matrix, which is stored in self.itof, gives
            %the conductivity in the forward mesh basis by multiplication
            %self.itof*sigma_i, where sigma_i is the conductivity in the
            %inverse mesh basis.
            if size(imesh.g,2) == 2%Check the dimension of the inverse mesh
                self.itof = self.interpolatematrix2d(imesh.H, imesh.g, self.g(:,1:2));
            elseif size(imesh.g,2) == 3
                self.itof = self.interpolatematrix3d(imesh.H, imesh.g, self.g);
            else
                error(['Unfit second dimension of ginv: size(ginv,2) = ' num2str(size(ginv,2))]);
            end
            self.nginv = size(imesh.g,1);
            
        end
        
        function SetEC(self)
            %Sets self.EC, i.e. an array, for which each row i contains all
            %such element indices (i.e. rows of H) that the element contains the node i (i.e. row i of g).
            %currently the only use for self.EC is by
            %self.GradientMatrix(). self.EC contains zeros on some of the
            %rows, as the number of elements touching a node is different
            %for some nodes, but self.EC is a matrix, so each row has the
            %same length

            counter = zeros(self.ng, 1);%a counter for how many elements are already mapped for each node
            self.EC = zeros(self.ng, 4);
            for iH = 1:self.nH%go through all elements
                for in = 1:size(self.H,2)%go through all the nodes of element iH
                    id = self.H(iH, in);
                    counter(id) = counter(id) + 1;%how many elements are already mapped to node id
                    self.EC(id, counter(id)) = iH;%map element iH to node id.
                end
            end
            
        end
        
        function int=intLinSigma(self,g,s,ip,L)
            %This function computes integrals of
            %s*grad(phi_i) dot grad(phi_j), where phi_i and phi_j are the
            %basis functions of triangular element vertices with
            %coordinates given by g. ip and L are auxiliary variables given
            %so that they do not have to be computed again for each element.

            %compute the conductivity values at the Gauss quadrature
            %integration points:
            sigma = [1-sum(ip,2) ip]*s;
            

            Jt = L*g;
            iJt = inv(Jt);
            dJt = abs(det(Jt));
            G = iJt*L;
            GdJt = G'*G*dJt/factorial(self.gdim+1);
            int = sum(sigma)*GdJt;
            %int is a 3-by-3 or 4-by-4 matrix whose element (i,j) is the
            %integral of s*grad(phi_i) dot grad(phi_j).

        end

    end

    methods(Static)

        function [int1, int2] = phiiphij1D(g)
            %This function computes the 1D integrals of phi_i and
            %phi_i*phi_j on an element defined by the coordinates in g
            
            a = 1/6;
            b = 1/3;
            c = 1/2;
            l=sqrt((g(1,1)-g(2,1))^2 + (g(1,2)-g(2,2))^2);
            int1 = l*[c;c];
            int2 = l*[b a;a b]; 

        end

        function [int1, int2] =phiiphij2D(g)
            %This function computes the 2D integrals of phi_i and
            %phi_i*phi_j on an element defined by the coordinates in g
            
            a = 1/24;
            b = 1/12;
            c = 1/6;
            A= norm(cross(g(2,:)-g(1,:), g(3,:)-g(1,:)));
            int1 = A*[c;c;c];
            int2 = A*[b a a;a b a; a a b]; 

        end

        function A=ElementArea2D(g)
            %Compute the are of element defined by the coordinates in g
            A= 0.5*norm(cross(g(2,:)-g(1,:), g(3,:)-g(1,:)));
        end

        function A=ElementArea1D(g)
            %Compute the are of element defined by the coordinates in g
            A = sqrt((g(1,1)-g(2,1))^2 + (g(1,2)-g(2,2))^2);
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