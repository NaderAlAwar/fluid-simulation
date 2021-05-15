import sys

assert(len(sys.argv) == 5)
nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])
mode = sys.argv[4] # stairs_1, stairs_2, fountain

def generateBarriersStairs1(nx: int, ny: int, nz: int) -> None:
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if i + j < ny // 2 and i + j >= (ny // 2 ) - 1:
                    print(f"c {i} {j} {k}")


def generateBarriersStairs2(nx: int, ny: int, nz: int) -> None:
    x, y = nx // 2, ny // 2
    count: int = 0
    while y >= 0 and x >= 0:
        if count < 4:
            for z in range(nz):
                print(f"c {x} {y} {z}")
                print(f"c {nx - x - 1} {y} {z}")
            x -= 1
        else:
            count %= 5
            y -= 1

        count += 1

def generateBarriersFountain(nx: int, ny: int, nz: int) -> None:
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if i + j < ny // 2 and i + j >= (ny // 2 ) - 1:
                    print(f"c {i} {j} {k}")
                    print(f"c {nx - i - 1} {j} {k}")

if mode == "stairs_1":
    generateBarriersStairs1(nx, ny, nz)
elif mode == "stairs_2":
    generateBarriersStairs2(nx, ny, nz)
elif mode == "fountain":
    generateBarriersFountain(nx, ny, nz)