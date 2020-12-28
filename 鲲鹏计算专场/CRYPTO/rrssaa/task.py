from Crypto.Util.number import *
from secret import flag 
from gmpy2 import next_prime,invert,gcd
import os                              

e = 65537
## level1
x = getPrime(290)
y = next_prime(21*x)
z = next_prime(3*x*y)
n1 = x*y*z
msg = flag+os.urandom(100)
m = bytes_to_long(msg)
assert(m < n1)
print n1
c1 = pow(m,e,n1)

## level2
m = c1
o = getPrime(300) 
s = getPrime(300)
t = next_prime(o)
u = next_prime(s)
print o*s
n2 = o*s*t*u
assert(m<n2)
print n2
c2 = pow(m,e,n2)
## level3
m = c2
p = getPrime(800)
q = getPrime(800)

n3 = p * q
phi = (q-1)*(p-1)
assert(m < n3)
print n3
while True:
    s = getPrime(10)
    if(gcd(s,p-1) == 1):
        sinv = invert(s,p-1)
        e = 4*s*sinv+3
        if(gcd(phi,e) == 1):
            break
c3 = pow(m,e,n3)
print c3
m = bytes_to_long(os.urandom(128))

print m
assert(m<n3)
print pow(m,e,n3)
