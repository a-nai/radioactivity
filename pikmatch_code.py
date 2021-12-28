def pikmatch(x,y,maxx,maxy,optimal,minmaxx,miny,minx):
    uj=0;pikk=[[float("inf") for m in range(len(maxy))] for t in range(len(maxx))];pi=[0 for m in range(len(maxx))];ret=1;p=[];nad=0;wax=[];way=[];re=1;
    for e in range(len(maxx)):
     for l in range(len(maxy)):  # b,n
      pikk[e][l]=(((abs(x[maxx[e]]-y[maxy[l]])+0.001)**1)*(abs(maxx[e]-maxy[l]+optimal)+1));
     pi[e]=min(pikk[e]);p.append(pikk[e].index(min(pikk[e])));p[e]=maxy[p[e]];nad=nad+pi[e];
    if len(minx)==0: minx=[];minx.append(1)
    if len(miny)==0: miny=[];miny.append(1)
    dob=1;res=0;pikk1=[[float("inf") for m in range(len(miny))] for t in range(len(minx))];p1=[];pi1=[0 for m in range(len(minx))];
    for e in range(len(minx)):
     for l in range(len(miny)):  # b,n
      pikk1[e][l]=(((abs(x[minx[e]]-y[miny[l]])+0.001)**1)*(abs(minx[e]-miny[l]-optimal)+1));
     pi1[e]=min(pikk1[e]);p1.append(pikk1[e].index(min(pikk1[e])));p1[e]=miny[p1[e]];nad=nad+pi1[e];
    dt1=dtw(pikk);dt=dt1[2];dt2=dtw(pikk1);dt3=dt2[2];
    for v in range(len(dt[1])): 
     p[dt[0][v]]=maxy[dt[1][v]];pi[dt[0][v]]=pikk[dt[0][v]][dt[1][v]];
     if (v>0)and(dt[0][v]==dt[0][v-1])and(pikk[dt[0][v-1]][dt[1][v-1]]<pi[dt[0][v]]): p[dt[0][v]]=maxy[dt[1][v-1]];pi[dt[0][v]]=pikk[dt[0][v-1]][dt[1][v-1]];
    for v in range(len(dt3[1])): 
     p1[dt3[0][v]]=miny[dt3[1][v]];pi1[dt3[0][v]]=pikk1[dt3[0][v]][dt3[1][v]];
     if (v>0)and(dt3[0][v]==dt3[0][v-1])and(pikk1[dt3[0][v-1]][dt3[1][v-1]]<pi1[dt3[0][v]]): p1[dt3[0][v]]=miny[dt3[1][v-1]];pi1[dt3[0][v]]=pikk1[dt3[0][v-1]][dt3[1][v-1]];
    nm=p[:];nm1=p1[:];rmu=1;pp=p;pp1=p1;
    D = defaultdict(list);D1 = defaultdict(list);
    for i,item in enumerate(p): D[item].append(i)
    for i,item in enumerate(p1): D1[item].append(i)
    D = {k:v for k,v in D.items() if len(v)>1};D1 = {k:v for k,v in D1.items() if len(v)>1} 
    #print p,p1,7#,minmaxx,optimal,pikk
    for v in range(len(list(D))):
     mah=inf;
     for e in range(len(D[list(D)[v]])): mah=pi[D[list(D)[v]][e]] if pi[D[list(D)[v]][e]]<mah else mah 
     k=pi.index(mah);
     if (minmaxx.index(maxx[k])+1)<len(minmaxx): perem=minx.index(minmaxx[minmaxx.index(maxx[k])+1])
     if (minmaxx.index(maxx[k])-1)>-1: perem1=minx.index(minmaxx[minmaxx.index(maxx[k])-1])
     for e in range(len(D[list(D)[v]])):
      if D[list(D)[v]][e]!=k:
       ms=D[list(D)[v]][e];p[ms]=p[k]-(maxx[k]-maxx[ms]);
       if ((maxx[k]-maxx[ms])<0)and(pp1.count(pp1[perem])==1): p1[perem]=p[ms]-(maxx[ms]-minx[perem]);
       if ((maxx[k]-maxx[ms])>0)and(pp1.count(pp1[perem1])==1): p1[perem1]=p[ms]-(maxx[ms]-minx[perem1]); 
    for v in range(len(D1.keys())):
     mah=inf;
     for e in range(len(D1[list(D1)[v]])): mah=pi1[D1[list(D1)[v]][e]] if pi1[D1[list(D1)[v]][e]]<mah else mah 
     k=pi1.index(mah);#print minmaxx,minx,maxx,(minmaxx.index(minx[k])+1),maxx.index(minmaxx[minmaxx.index(minx[k])+1])
     if (minmaxx.index(minx[k])+1)<len(minmaxx): perem=maxx.index(minmaxx[minmaxx.index(minx[k])+1])
     if (minmaxx.index(minx[k])-1)>-1: perem1=maxx.index(minmaxx[minmaxx.index(minx[k])-1])
     for e in range(len(D1[list(D1)[v]])):
      if D1[list(D1)[v]][e]!=k:
       ms=D1[list(D1)[v]][e];p1[ms]=p1[k]-(minx[k]-minx[ms]);
       if ((minx[k]-minx[ms])<0)and(pp.count(pp[perem])==1): p[perem]=p[ms]-(minx[ms]-maxx[perem]);
       if ((minx[k]-minx[ms])>0)and(pp.count(pp[perem1])==1): p[perem1]=p[ms]-(minx[ms]-maxx[perem1]);
    for i in range(len(p1)): 
     if p1[i]>59: p1[i]=59; 
     if p1[i]<0: p1[i]=0; 
    for i in range(len(p)): 
     if p[i]>59: p[i]=59; 
     if p[i]<0: p[i]=0; 
    #print p,p1,5
    for e in range(len(p)+len(p1)):  #+x[maxx[l]]+y[p[e]]+y[p[l]]     ,len(p)+len(p1)
     for l in range(len(p)+len(p1)):  #-e +len(p1)
      if e>=len(p): me=minx[e-len(p)];pe=p1[e-len(p)];nme=nm1[e-len(p)];
      else: me=maxx[e];pe=p[e];nme=nm[e];
      if l>=len(p): ml=minx[l-len(p)];pl=p1[l-len(p)];nml=nm1[l-len(p)];
      else: ml=maxx[l];pl=p[l];nml=nm[l];
      red1=1;red=1;ab=y[pl]-y[pe];ac=pe-pl;au=x[ml]-x[me];ao=me-ml;hip=math.sqrt(ab**2+ac**2);hipe=math.sqrt(ao**2+au**2);
      if (hip!=0)and(hipe!=0): 
       if math.asin(abs(au)/hipe)!=0:
        red=abs(abs(math.asin(abs(ab)/hip)/math.asin(abs(au)/hipe)));
      #if red==0: red=1;
      #if red<1: red=1.0/red;#print red;
      ret=abs(1-abs((0.00+abs((y[pl]+0.01)/(y[pe]+0.01)))/(0.00+abs((x[ml]+0.01)/(x[me]+0.01)))));uj=uj+1;
      re=0+abs(1-abs(float(abs(pe-pl)+1))/float(abs(me-ml)+1));rme=(y[pe]+1)*(y[pl]+1)*(x[me]+1)*(x[ml]+1);rmu+=rme; #(x[me]+1)*(x[ml]+1)*
      res+=abs(re*rme*red*ret);
    return len(pi)*len(pi),res/uj; 
