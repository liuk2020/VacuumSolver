from vacuum.fortraninterface import get_resolution


xm, xn = get_resolution(mpol=2, ntor=3, nfp=5)
print(xm)
print(xn)