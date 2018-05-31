function circle(r,side)
% draws a filled circle with fill defined by side
% r: defines radius
% side: defines the side (1= fill outside, 0=foll inside)
    t=linspace(0,2*pi,1000);
    x=r*cos(t); y=r*sin(t);   % create circle points for ROC
    if side==1
        c=[1 1 1];
        set(gca,'color',[.95 .95 .95]);
        hold(gca,'on');
    elseif side==2 % fill circle white and white axes
        c=[1 1 1];
    else
        c=[.95 .95 .95];
    end
    fill(x,y,c); % draw circle    
    hold(gca,'off');
end