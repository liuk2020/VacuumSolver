from vacuum.fortraninterface import set_resolution


xm, xn = set_resolution(mpol=2, ntor=3, nfp=5)
print(xm)
print(xn)