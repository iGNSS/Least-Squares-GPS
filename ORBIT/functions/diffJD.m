function s = diffJD(t1, t2)

%difference t1 - t2
dd = t1.d - t2.d;
ds = t1.s - t2.s;
s = dd*86400 + ds;
