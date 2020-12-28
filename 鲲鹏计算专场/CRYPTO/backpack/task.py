import random
import hashlib
from functools import reduce
from binascii import hexlify, unhexlify
from Crypto.Util.number import getPrime, inverse, long_to_bytes, bytes_to_long
from secret import flag, nonce1, nonce2

def crt(n, a):
    sum = 0
    prod = reduce(lambda a, b: a * b, n)
    for n_i, a_i in zip(n, a):
        p = prod // n_i
        sum += a_i * inverse(p, n_i) * p
    return int(sum % prod)

class Encryptor:
    def __init__(self, n):
        self.u = []
        self.v = []
        self.pubkey = []
        self.n = n

    def generate(self, ):
        self.u = [random.randint(1, 64) for i in range(self.n)]
        self.v = [self.u[i] - pow(2, self.n - i - 1) for i in range(self.n)]
        t1 = sum(self.u)
        t2 = 2 * max(sum(self.v), -sum(self.v))
        while 1:
            self.p = getPrime(t1.bit_length() + 10)
            self.q = getPrime(t2.bit_length() + 10)
            if self.p > t1 and self.q > t2:
                self.N = self.p * self.q
                for i in range(self.n):
                    self.pubkey.append(crt([self.p, self.q], [self.u[i], self.v[i]]))
                break

    def encrypt(self, msg):
        tmp = bin(bytes_to_long(msg))[2:].zfill(self.n)
        c = 0
        for i in range(self.n):
            c += int(tmp[i]) * self.pubkey[i]
        return c


assert hashlib.sha1(nonce1).hexdigest() == "51d6169bcc32acb2a4d3b1a8d9c6ed0c9a909974"
assert hashlib.sha1(nonce2).hexdigest() == "2347411264fc395375fdfe3dbd6169283f3e4923"

encryptor1 = Encryptor(100)
encryptor1.generate()
ct1 = encryptor1.encrypt(nonce1)
print(f"pubkey:{encryptor1.pubkey}")
print(f"ct1:{ct1}")

encryptor2 = Encryptor(320)
encryptor2.generate()
ct2 = encryptor2.encrypt(nonce2)
print(f"pubkey:{encryptor2.pubkey}")
print(f"ct2:{ct2}")

assert flag == "flag{" + hashlib.sha256(nonce1).hexdigest()[:16] + hashlib.md5(nonce2).hexdigest()[:16] + "}"

