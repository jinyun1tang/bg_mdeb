function [t,y]=euler_driver_sdeb(y0,tspan,ns)

dt=(tspan(2)-tspan(1))/ns;
y1=y0;
t=(tspan(1):dt:tspan(2));
nel=numel(t);
n=0;
y=zeros(nel,numel(y0));
nl=numel(y0);
while(n<nel)
    y(n+1,:)=y1;
    dydt2=sDEB(t(n+1),y1);
    y2=y1+dydt2.*dt;
    pp=1;
    for jj = 1 : nl
        if(y2(jj)<0)
            pp1=abs(y1(jj)/(dydt2(jj)*dt))*0.9999;
            pp=min(pp,pp1);
        end
    end
    if(pp<1)
        y1=y1+dydt2.*dt.*pp;
    else
        y1=y2;
    end
    n=n+1;
end
end