function p = curvy(coefs,subd) 

%input control points, returns the nurbs curve points
%
%Default value of knots = [0,0,0,0,1,1,1,1]
%Default value of subd = 100
%change the default values by editing this function

knots = [0,0,0,0,1,1,1,1];
% subd = 100;

nurbs.form   = 'B-NURBS'; 
nurbs.dim    = 4; 
np = size(coefs); 
dim = np(1); 

% constructing a curve 
nurbs.number = np(2); 
if (dim < 4) 
  nurbs.coefs = repmat([0.0 0.0 0.0 1.0]',[1 np(2)]); 
  nurbs.coefs(1:dim,:) = coefs;   
else 
  nurbs.coefs = coefs; 
end 
nurbs.order = size(knots,2)-np(2); 
knots = sort(knots); 
nurbs.knots = (knots-knots(1))/(knots(end)-knots(1));

tt = linspace(0.0,1.0,subd);

d = nurbs.order-1;
c = nurbs.coefs;
k = nurbs.knots;
u = tt;


nu = numel(u); 
[mc,nc] = size(c); 
                                                                                                                                                                                                                                                                                                                                         
p = zeros(mc,nu);                                
N = zeros(d+1,1);                               
                                                 
                                               
for col=1:nu                                   
                                                
  s = findspan(nc-1, d, u(col), k);            
  N = basisfun(s,u(col),d,k);                 
                                                
  tmp1 = s - d + 1;                           
  for row=1:mc                                
      tmp2 = 0;                               
      for i=0:d                               
         tmp2 = tmp2 + N(i+1)*c(row,tmp1+i);   
      end                                      
      p(row,col) = tmp2;                      
  end
end


w = p(4,:); 
p = p(1:3,:); 
p = p./repmat(w,3,1); 

end 
 