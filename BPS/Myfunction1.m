function fun = Myfunction1( )

a=1;
a2=1;
a3=0.001;


fun{1,1} = @(x)  (x.^2-a).^2/4 ;
fun{1,2} = @(x)  x.^4/4 ;

% Energy
fun{1,3} = @(x) (x(1).^2-a ).^2/4+ x(2).^4/4 ;
fun{1,4} = @(x,y) (x.^2-a ).^2/4+ y.^4/4 ;
%  Force
fun{2,1} = @(x) a3*[x(1).^3-a2*x(1);x(2).^3];
% Nabla Force
fun{2,2} = @(x)  a3*[3*x(1).^2-a2,0;0,3*x(2).^2];

%  nabla U
fun{2,3} = @(x) [x(1).^3-a*x(1);x(2).^3];

% nabla nabla U
fun{2,4} = @(x)  [3*x(1).^2-a,0;0,3*x(2).^2];
