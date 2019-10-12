function [] = plot_theta(theta,pp)

%% Draw circle
style='k'; NOP=100; 
center = [0 0]; radius = 100;
THETA=linspace(0,(2*pi),NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
plot(X,Y,style);
axis square;

%% Plot angles
xx = radius*cos(theta); 
yy = radius*sin(theta);
hold on; scatter(xx,yy,'b');

if ~isempty(pp)
    xxp=radius*cos(pp);
    yyp=radius*sin(pp);
    scatter(xxp,yyp,'r.')
    
    xxp=radius*cos((pi/2)+(pi/2)-pp);
    yyp=radius*sin((pi/2)+(pi/2)-pp);
    scatter(xxp,yyp,'r.')
    
    xxp=radius*cos(0.5*(pp+((pi/2)+(pi/2)-pp)));
    yyp=radius*sin(0.5*(pp+((pi/2)+(pi/2)-pp)));
    scatter(xxp,yyp,'ks')
    
end


plot(0,linspace(-100,100,200),'k-'); plot(linspace(-100,100,200),0,'k-')
text(-50,60,'D','Color','black','FontSize',16,'FontName','Arial','fontweight','bold')
text(50,60,'T','Color','black','FontSize',16,'FontName','Arial','fontweight','bold')

hold off;
end