import subprocess
import hashlib
import base64
import os

from shutil import copyfile

flag = "flag{test}"

def check_output():
    f1 = open("/root/check.dat", "rb")
    f2 = open("/tmp/data.dat", "rb")
    data1 = np.fromfile(f1, dtype=np.float32)
    data2 = np.fromfile(f2, dtype=np.float32)
    if  ((data1-data2) < 0.00001).all():
        print(flag)

if __name__ == "__main__":
    code = input("Encode your executable using base64: ")
    print()
    if not code:
        print("Code must not be empty")
        exit(-1)
    elf = base64.b64decode(code)
    with open("/tmp/elf","wb") as f:
        f.write(elf)
        
    os.chmod("/tmp/elf", 0o755)
    copyfile("/root/DATA.bin", "/tmp/DATA.bin")
    p = subprocess.run(
        ["su", "nobody", "-s", "/tmp/elf"],
        input=elf_input.encode(),
        stdout=subprocess.PIPE,
    )

    if p.returncode != 0:
        print()
        print("Your code did not run successfully")
        exit(-1)
    check_output()