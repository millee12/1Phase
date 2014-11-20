    figure;
    for a=1:nf
        i=2;
        X=0;
        Y=0;
        X(1)=state.mesh(a);
        for b=1:ns
            if Ind(a,b) > 0
                if sum(Ind(a,1:b)) < 1
                    X(i)=state.geometry.breakpoints(b+1);
                    Y(i)=state.geometry.crosssections(b);
                    X(i+1)=state.geometry.breakpoints(b+1);
                    Y(i+1)=state.geometry.crosssections(b+1);
                    i=i+2;
                else
                    X(i)=state.mesh(a+1);
                    Y(i)=state.geometry.crosssections(b);
                    i=i+1;
                end
            end
        end
        Y(1)=Y(2);
        X=[X,X(end:-1:1)];
        Y=[Y,-Y(end:-1:1)];
        figure(gcf);
        p1=patch(X,Y,[0,0,0]);
        set(p1,'FaceColor','none','EdgeColor','b','Linewidth',2)
    end