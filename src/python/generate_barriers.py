import sys

assert(len(sys.argv) == 4)
nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])

def generateBarriers(nx: int, ny: int, nz: int) -> None:
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if i + j < ny // 2 and i + j >= (ny // 2 ) - 1:
                    print(f"c {i} {j} {k}")

generateBarriers(nx, ny, nz)