def stat_distanc(counter):
    koeff=60;ex=20;ko=60;b=[0 for e in range(ko+ex)];bb=[0 for e in range(ko+ex)];#C:\\28-08-2004.txt /usr/local/28-08-2004.txt
    matrix = [line.strip() for line in open('/users/andrejeremcuk/downloads/0505.txt')]; #days.dat
    matri = [line.strip() for line in open('/users/andrejeremcuk/downloads/0506.txt')]; #days.dat
    arra=[[ 0 for e in range(koeff)] for t in range(int(len(matrix)/koeff))];harra=[0 for t in range(int(len(matrix)/koeff))];
    arrab=[[ 0 for e in range(koeff)] for t in range(int(len(matrix)/koeff))];harrab=[0 for t in range(int(len(matrix)/koeff))];z=0;resume=1;resum=0;
    arr=[[ 0 for e in range(koeff)] for t in range(int(len(matri)/koeff))];
    for i in range(len(matrix)): matrix[i]=int(matrix[i]);matri[i]=int(matri[i]);
    for jk in range(2880): #int(len(matrix)/koeff)
     for mk in range(koeff): arra[jk][mk]=matrix[z];arr[jk][mk]=matri[z];z=z+1;
    for jk in range(2880):
     harra[jk]=np.histogram(arra[jk], bins=ko, density=False);harrab[jk]=np.histogram(arr[jk], bins=ko, density=False);harra[jk]=harra[jk][0];harrab[jk]=harrab[jk][0]
     for u in range(5): harra[jk] = movingaverage(harra[jk], 4);harrab[jk]=movingaverage(harrab[jk], 4);
    for mk in range(1440):# 24000
     mi=[];mis=[];mis1=[];mii=[];ei=[];ti=[];dots=[];dots1=[];
     if counter=='r': rn=random.randint(0,1439);x=harra[mk];y=harrab[random.randint(0,1439)];#print rn;  #y=harra[random.randint(0,1439)];
     elif counter=='o': x=harra[summ[mk]-1];y=harra[sum1[mk]-1];
     else: x=harra[mk+1440-counter];y=harrab[mk];
     x=x.tolist();y=y.tolist();
     for k in range(ko): b[int(k+ex/2)]=y[k];
     for m in range(ex):
      ct=[0 for e in range(ko+ex)]; 
      for l in range(ko): ct[l+m]=x[l];
      mis.append(euclid(ct,b));mis1.append(euclid(ct[::-1],b));
     optimal=mis.index(min(mis))-ex/2;ret=pik(x,y);minmaxx=ret[6];minmaxy=ret[7];maxx=ret[0];maxy=ret[1];minx=ret[2];miny=ret[3]
     optimal1=mis1.index(min(mis1))-ex/2;#print mk;
     for k in range(len(minmaxx)): dots.append(x[minmaxx[k]]);#minmaxx[k]=59-minmaxx[k];
     for k in range(len(minmaxy)): dots1.append(y[minmaxy[k]]);#minmaxy[k]=59-minmaxy[k]
     outputs=pikmatch(x,y,maxx,maxy,-optimal,minmaxx,miny,minx);
     outputs1=pikmatch(y,x,maxy,maxx,optimal,minmaxy,minx,miny);
     #.set_alpha(0.5)
     mx1=minmaxx[:];my1=minmaxy[:];mxx1=maxx[:];mxx2=minx[:]
     for k in range(len(minmaxx)): dots.append(x[minmaxx[k]]);minmaxx[k]=59-minmaxx[k];
     for k in range(len(maxx)): maxx[k]=59-maxx[k];
     for k in range(len(minx)): minx[k]=59-minx[k];
     minmaxx=minmaxx[::-1];maxx=maxx[::-1];minx=minx[::-1];
     output=pikmatch(x[::-1],y,maxx,maxy,-optimal1,minmaxx,miny,minx);minmaxx=minmaxx[::-1];mx2=minmaxx[:]; #
     for k in range(len(minmaxy)): dots1.append(y[minmaxy[k]]);#minmaxy[k]=59-minmaxy[k]
     minmaxx=mx1[:];
     output1=pikmatch(y,x[::-1],maxy,maxx,optimal1,minmaxy,minx,miny);
     resume+=min(outputs1[1]+outputs[1],output[1]+output1[1]);
     #minmaxx=mx2[:];plt.plot(minmaxx,dots,'go');plt.plot(my1,dots1,'bo');plt.plot(y);plt.plot(x[::-1]);plt.show(); #plt.plot(x);
    return resume/10000;

stat_distanc('r')
stat_distanc(0)

stat_distanc(1072)

bas=[];
for i in range(1440): bas.append(stat_distanc(i));print(i);
    
