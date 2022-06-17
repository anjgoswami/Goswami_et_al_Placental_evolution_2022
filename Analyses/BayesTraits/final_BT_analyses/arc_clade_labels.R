require(plotrix)
arc.cladelabels<-function(tree=NULL,text,node,ln.offset=1.02,
    lab.offset=1.06,cex=1,orientation="curved",...){
    obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
    if(obj$type!="fan") stop("method works only for type=\"fan\"")
    h<-max(sqrt(obj$xx^2+obj$yy^2))
    if(hasArg(mark.node)) mark.node<-list(...)$mark.node
    else mark.node<-TRUE
    if(mark.node) points(obj$xx[node],obj$yy[node],pch=21,
        bg="red")
    if(is.null(tree)){
        tree<-list(edge=obj$edge,tip.label=1:obj$Ntip,
            Nnode=obj$Nnode)
        class(tree)<-"phylo"
    }
    d<-getDescendants(tree,node)
    d<-sort(d[d<=Ntip(tree)])
    deg<-atan(obj$yy[d]/obj$xx[d])*180/pi
    ii<-intersect(which(obj$yy[d]>=0),which(obj$xx[d]<0))
    deg[ii]<-180+deg[ii]
    ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]<0))
    deg[ii]<-180+deg[ii]
    ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]>=0))
    deg[ii]<-360+deg[ii]
    draw.arc(x=0,y=0,radius=ln.offset*h,deg1=min(deg),
        deg2=max(deg))
    if(orientation=="curved")
        arctext(text,radius=lab.offset*h,
            middle=mean(range(deg*pi/180)),cex=cex)
    else if(orientation=="horizontal"){
        x0<-lab.offset*cos(median(deg)*pi/180)*h
        y0<-lab.offset*sin(median(deg)*pi/180)*h
        text(x=x0,y=y0,label=text,
        adj=c(if(x0>=0) 0 else 1,if(y0>=0) 0 else 1),
        offset=0)
    }
}