import subprocess
import hashlib
import base64
import os

import string
import random
from Crypto.Cipher import AES

flag = "flag{test}"

key = ''.join(random.choices(string.ascii_letters+string.digits,k=16))

def gen_input():
    cipher = AES.new(key.encode(), AES.MODE_ECB)
    ct = cipher.encrypt(b'KUNPENG_HPC_AES!')
    return key[:12]+ct.hex()

def check_output(elf_input, output):
    if elf_input[:12] + output == key:
        return True
    return False

if __name__ == "__main__":
    code = input("Encode your executable using base64: ")
    print()
    if not code:
        print("Code must not be empty")
        exit(-1)
    elf = base64.b64decode(code)
    with open("/elf","wb") as f:
        f.write(elf)
        
    os.chmod("/elf", 0o755)
    elf_input = gen_input()
    print("the input is: " + elf_input)
    p = subprocess.run(
        ["su", "nobody", "-s", "/elf"],
        input=elf_input.encode(),
        stdout=subprocess.PIPE,
    )

    if p.returncode != 0:
        print()
        print("Your code did not run successfully")
        exit(-1)

    output = p.stdout.decode()
    print("Your output is: " + output)
    if check_output(elf_input, output):
        print(flag)
    else:
        print("Wrong")