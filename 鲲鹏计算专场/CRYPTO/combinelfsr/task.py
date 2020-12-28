from secret import flag
from hashlib import sha512
import random
from Crypto.Cipher import AES

def lfsr(R,mask):
        output = (R << 1) & 0xffffff
        i=(R&mask)&0xffffff
        lastbit=0
        while i!=0:
            lastbit^=(i&1)
            i=i>>1
        output^=lastbit
        return (output,lastbit)
up = 2**18
R1 = random.randint(0,up)
R2 = random.randint(0,up)
mask1 = 0x30517
mask2 = 0x25b74
print R1
print R2

key = sha512(str(R1)+str(R2)).digest()[:16]
aes = AES.new(key,AES.MODE_ECB)
flag += '\x00'*(16-len(flag)%16)
flag = aes.encrypt(flag).encode('hex')
f = open('flag.txt','wb')
f.write(flag)
f.close()

f = open('out','wb')
def combine(r1,r2,mask1,mask2):
    (r11,x1)=lfsr(r1,mask1)
    (r22,x2)=lfsr(r2,mask2)
    return (r11,r22,(x1*x2)^(x2^1))
for i in range(1024):
    tmp=0
    for k in range(8):
        R1,R2,out = combine(R1,R2,mask1,mask2)
        tmp = (tmp << 1) ^ out
    f.write(chr(tmp))
f.close()



