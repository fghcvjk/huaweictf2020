import subprocess
import hashlib
import base64
import os

import string
import random

flag = "flag{test}"

def gen_input():
    s=string.ascii_letters+string.digits
    return ''.join(random.choices(s,k=10))

def check_output(elf_input, output):
    if bin(int.from_bytes(hashlib.sha256((elf_input+output).encode()).digest(), "big"))[2:].zfill(256).startswith('0'*30):
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