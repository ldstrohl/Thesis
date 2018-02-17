%---------------------------------------------------------------
function [za]=eqn535(xa,ya)
% limit inputs
if xa > 3
    xa = 3;
elseif xa < 0.2
    xa = 0.2;
end

if ya > 90
    ya = 90;
elseif ya < 0
    ya = 0;
end

%---------------------------------------------------------------
%   TableCurve 3D
%   File Source= d:\sdsu\graduate students and postdoc\lioyd strohl\mrv\cl_supersonic-v2.dat
%   Date= Aug 9, 2017
%   Time= 9:18:38 PM
%   Data Set= CL_supersonic-v2.dat, X , Y , Z 
%   X= 
%   Y= 
%   Z= 
%   Eqn#= 535
%   Eqn= Cosine Series Bivariate Order 6
%   r2=0.9751427530542483
%   r2adj=0.9702413240790297
%   StdErr=0.138824367276737
%   Fstat=207.7722027154769
%   a= 1.209533073301437
%   b= 0.1615343388239877
%   c= -0.126993522475366
%   d= -0.1258635264316632
%   e= -0.005720656579217343
%   f= -0.9880503149453188
%   g= -0.1013594115014437
%   h= -0.08389605072166674
%   i= -0.02994434949806103
%   j= 0.135108330843908
%   k= -0.03452450639711589
%   l= -0.05298715003955762
%   m= 0.1194976763053231
%   n= 0.01827696843647851
%   o= -0.1455900663246896
%   p= -0.01530945350276446
%   q= -0.02781412939505948
%   r= 0.1480066199295564
%   s= 0.04746803879093773
%   t= -0.01347573690331952
%   u= 0.02091090804661737
%   v= 0.03093510610323948
%   aa= -0.02777429429730163
%   ab= 0.1617055307517038
%   ac= 0.03171620146620202
%   ad= -0.04165982604234347
%   ae= -0.02236515087593543
%   af= -0.07580602962244646
  [rowx colx]=size(xa);
  if(rowx~=1 & colx~=1)
    error('x must be scalar or 1D array');
    return;
  end
  [rowy coly]=size(ya);
  if(rowy~=1 & coly~=1)
    error('y must be scalar or 1D array');
    return;
  end
  c=[
    1.209533073301437,
    0.1615343388239877,
    -0.1269935224753660,
    -0.1258635264316632,
    -0.005720656579217343,
    -0.9880503149453188,
    -0.1013594115014437,
    -0.08389605072166674,
    -0.02994434949806103,
    0.1351083308439080,
    -0.03452450639711589,
    -0.05298715003955762,
    0.1194976763053231,
    0.01827696843647851,
    -0.1455900663246896,
    -0.01530945350276446,
    -0.02781412939505948,
    0.1480066199295564,
    0.04746803879093773,
    -0.01347573690331952,
    0.02091090804661737,
    0.03093510610323948,
    -0.02777429429730163,
    0.1617055307517038,
    0.03171620146620202,
    -0.04165982604234347,
    -0.02236515087593543,
    -0.07580602962244646,
    ];
  lenx=length(xa);
  leny=length(ya);
  for(j=1:leny)
    for(i=1:lenx)
      x=xa(i);
      y=ya(j);
      z=evalcsi(27,x,y,c,...
        0.2000000000000000,0.8912676813146139,...
        0.000000000000000,28.64788975654116);
       za(i,j)=z;
       end
     end
%--------------------------------------------------------------
function z = evalcsi(order, x, y, p, s0, s1, s2, s3)
%--------------------------------------------------------------
  cx=[];
  cy=[];
  v=[];
  x=(x-s0)/s1;
  y=(y-s2)/s3;
  if(x<0.0)
    x=0.0;
  end
  if(x>3.14159265358979323846)
    x=3.14159265358979323846;
  end
  if(y<0.0)
    y=0.0;
  end
  if(y>3.14159265358979323846)
    y=3.14159265358979323846;
  end
  if (order==5)
    nc=2;
  elseif (order==9)
    nc=3;
  elseif (order==14)
    nc=4;
  elseif (order==20)
    nc=5;
  elseif (order==27)
    nc=6;
  elseif (order==35)
    nc=7;
  elseif (order==44)
    nc=8;
  elseif (order==54)
    nc=9;
  elseif (order==65)
    nc=10;
  else
    return;
  end
  cx(1)=1.0;
  cy(1)=1.0;
  for(j=1:1:nc)
    cx(j+1)=cos(j*x);
    cy(j+1)=cos(j*y);
  end
  iv=1;
  for(j=0:1:nc)
    for(m=j:-1:0)
      v(iv)=cx(m+1)*cy(j-m+1);
      iv=iv+1;
    end
  end
  z=0.0;
  for(j=1:1:order+1)
    z=z+p(j)*v(j);
  end
  return;
