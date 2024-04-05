function deltareco_pixgrid = interpolateRecoToPixGrid(deltareco,Mesh);
R = range(Mesh.g(:,1))/2;
pixwidth = 2*R/256;
pixcenter_x = [-R+pixwidth/2:pixwidth:R-pixwidth/2];
pixcenter_y = pixcenter_x;
[X,Y] = meshgrid(pixcenter_x,pixcenter_y);
pixcenters = [X(:), Y(:)];
deltareco_pixgrid = Interpolate2Newmesh2DNode(Mesh.g,Mesh.H,[],deltareco,pixcenters,[]);
deltareco_pixgrid = flipud(reshape(deltareco_pixgrid,256,256));
end

%Interpolates the electric conductivity (discretized using linear basis
%functions) to a new mesh.
function [f_newgrid,INTPMAT,Element] = Interpolate2Newmesh2DNode(g,H,~,f,pts,INTPMAT)

if(isempty(INTPMAT))

   invX = cell(size(H,1),1);
   for tin = 1:size(H,1);

        pp = H(tin,:);
        X = [g(pp,:),ones(3,1)];
        invX{tin} = inv(X);
 
   end
   

   np = size(pts,1);
   Ic = zeros(np,3);
   Iv = zeros(np,3);

TR = triangulation(H,g);
   Element = pointLocation(TR,pts);
   
   nans = zeros(np,1);
   for k = 1:np

      x = pts(k,1); y = pts(k,2);
      tin = Element(k);
      
      Phi = zeros(1,3);
      if ~isnan(tin)
      if tin       
        iXt = invX{tin};
        for gin = 1:3
          Phi(gin) = [x y 1]*iXt(:,gin);
        end
        Ic(k,:) = H(tin,:);
      else  
      end
            
      Iv(k,:) = Phi;  
      else
          Ic(k,:) = 1;
          nans(k) = 1;
          Iv(k,:) = 1;  
      end
   end
   INTPMAT = sparse(repmat([1:np]',1,3),Ic,Iv,np,size(f,1));
   INTPMAT(nans==1,:) = 0;
   
end

f_newgrid = INTPMAT*f;

end




