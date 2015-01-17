delta=[1 1/2 1/4 1/8 1/16 1/32];
for i=1:6
    i
    x=Heston_SZ2(100,delta(i),27.9,4)
end